import MolecularGraph: 
    MolGraph,
    SDFAtom, SDFBond, 
    SMILESAtom, SMILESBond,
    Metadata, Stereobond,
    Stereocenter,
    has_state,
    set_state!

function _atom_to_molgraph(a)
    Dict{String, Any}(
        "stereo" => :unspecified, # TODO: handle stereo chemistry
        "multiplicity" => 1,      # TODO: handle multiplicities
        "symbol" => a.element,
        "charge" => a.formal_charge,
        "mass"   => nothing,      # TODO: handle masses
        "coords" => a.r,
        "idx"    => a.idx,
        "number" => a.number
    )
end

function _bond_to_molgraph_edge(b, idx_to_molgraph_atom)
    (idx_to_molgraph_atom[b.a1], idx_to_molgraph_atom[b.a2])
end

function _bond_to_molgraph_attr(b)
    Dict{String, Any}(
        "notation" => 0,            # TODO: handle notation
        "stereo"   => :unspecified, # TODO: handle stereo chemistry
        "order"    => b.order <= BondOrder.Quadruple ? Int(b.order) : 1, # TODO: handle aromatic bonds correctly
        "idx"      => b.idx,
        "isordered"=> true          # needed attribute for MolGraph SDFBonds
    )
end

function _gprops_dict_to_vec(d) 
    # serializes the "metadata", "stereocenter" and "stereobond" 
    # entries in the properties of a AbstractAtomContainer.
    # these entries are types that only carry a dict called "mapping"
    dummy = (mapping=Dict(),)
    subdicts = getproperty.((get(d, :metadata, dummy), get(d, :stereocenter, dummy), get(d, :stereobond, dummy)), :mapping)
    temp1 = Dict(k => string(v) for (k,v) in subdicts[1])
    temp2 = [vcat(k,string(v)) for (k,v) in subdicts[2]]
    temp3 = [vcat(e.src, e.dst, bond) for (e,bond) in subdicts[3]]
    props = Any[temp1, temp2, temp3]                            
    #isempty(props[1]) && (props[1] = Dict{String, Any}())           # MolecularGraphs wants this to be an empty dict if metadata is empty...
    return collect(zip( ("metadata", "stereocenter", "stereobond"),
                        string.((Metadata, Stereocenter{Int}, Stereobond{Int})), 
                        props))

    
end

function Base.convert(
        ::Type{MolGraph{SDFAtom, SDFBond}},
        mol::AbstractAtomContainer{T}
    ) where {T<:Real}

    # create an intermediate dictionary to map from our data structures
    # to those of MolecularGraph
    molgraph_atoms = map(_atom_to_molgraph, atoms(mol))
    idx_to_molgraph_atom = Dict(
        a["idx"] => i for (i,a) in enumerate(molgraph_atoms)    #todo:use ordered dict instead.
    )

    molgraph_atom_to_num = sort(molgraph_atoms, by=(at-> at["number"])) #need this to be sorted because of iter order for "graph"

    molgraph_atom_to_idx = Dict(
        v => k for (k,v) in idx_to_molgraph_atom
    )


    d = Dict{String, Any}(
        "vproptype"     => "SDFAtom",
        "vprops"        => molgraph_atom_to_num,        #graph atom => atominfo
        "eproptype"     => "SDFBond",
        "graph"         => map(b -> _bond_to_molgraph_edge(b, idx_to_molgraph_atom), bonds(mol)),
        "eprops"        => map(_bond_to_molgraph_attr, bonds(mol)),
        "caches"        => Dict{Any, Any}(),
        "gprops"        => _gprops_dict_to_vec(mol.properties)      #Todo: add here molgraph_idx => ball_idx dict
    )

    mol_graph = MolGraph{Int, SDFAtom, SDFBond}(d)
    set_state!(mol_graph, :atom_idx, molgraph_atom_to_idx)

    return mol_graph

end

function _molgraph_to_atom((i, a)::Tuple{Int, SDFAtom}, T)
    (
        number  = i,
        element = getproperty(Elements, a.symbol),
        name    = "$(a.symbol)$(i)",
        atomtype = "",
        r = Vector3{T}(a.coords),
        v = zeros(Vector3{T}),
        F = zeros(Vector3{T}),
        formal_charge = a.charge,
        charge = zero(T),
        radius = zero(T),
        properties = Properties(
            :multiplicity => a.multiplicity,
            :mass         => a.mass
        )
    )
end

function _molgraph_to_atom((i, a)::Tuple{Int, SMILESAtom}, T)
    (
        number  = i,
        element = getproperty(Elements, a.symbol),
        name    = "$(a.symbol)$(i)",
        atomtype = "",
        r = zeros(Vector3{T}),
        v = zeros(Vector3{T}),
        F = zeros(Vector3{T}),
        formal_charge = a.charge,
        charge = zero(T),
        radius = zero(T),
        properties = Properties(
            :multiplicity => a.multiplicity,
            :mass         => a.mass,
            :stereo       => a.stereo,
            :is_aromatic  => a.isaromatic
        )
    )
end

function _molgraph_to_bond((i, (e, b))::Tuple{Int, Tuple{Any, SDFBond}}, mol)
    df = atoms_df(mol)

    (
        a1 = df[df.number .== e.src, :idx][1],
        a2 = df[df.number .== e.dst, :idx][1],
        order = b.order,
        properties = Properties(
            :notation => b.notation
            )
    )
end

function _molgraph_to_bond((i, (e, b))::Tuple{Int, Tuple{Any, SMILESBond}}, mol)
    df = atoms_df(mol)

    (
        a1 = df[df.number .== e.src, :idx][1],
        a2 = df[df.number .== e.dst, :idx][1],
        order = b.order,
        properties = Properties(
            :is_aromatic => b.isaromatic,
            :direction   => b.direction,
            :stereo      => b.stereo
        )
    )
end

function Base.convert(
    ::Type{Molecule{T}},
    mg::MolGraph{GMAtom, GMBond};
    system=default_system()) where {T<:Real, GMAtom, GMBond}
    
     
    mol = Molecule(system, "", Dict{Symbol, Any}(Symbol(k) => v for (k, v) in mg.gprops))

    for a in (t -> _molgraph_to_atom(t, T)).(zip(keys(mg.vprops), values(mg.vprops)))
        Atom(mol, a...)
    end

    if has_state(mg, :atom_idx)
        for b in _molgraph_to_bond.(enumerate(( (src = mg.state[k.src], dst = mg.state[k.dst] ) ,v) for (k,v) in mg.eprops), Ref(mol))  
            Bond(parent_system(mol), b.a1, b.a2, BondOrderType(b.order), b.properties)
        end
    else
        for b in _molgraph_to_bond.(enumerate((k,v) for (k,v) in mg.eprops), Ref(mol))
            Bond(parent_system(mol), b.a1, b.a2, BondOrderType(b.order), b.properties)
        end
    end
    
    mol
end