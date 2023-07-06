using AutoHashEquals
using CSV
using Pipe
export MMFF94Parameters, _read_types_paramfile, assign_atom_types_mmff94,
assign_bond_types_mmff94, assign_charges

###### Constants #######
bend_elems = (Elements.H, Elements.B,Elements.C,Elements.N,Elements.O,Elements.F,
    Elements.Si,Elements.P,Elements.S,Elements.Cl,Elements.Br,Elements.As,Elements.I)
bend_z_ = (1.395, 0. , 2.494, 2.711, 3.045, 2.847, 2.350, 2.350, 2.980, 2.909, 3.017, 0., 3.086)
bend_c_ = (0.,    0.704, 1.016, 1.113, 1.337, 0., 0.811, 1.068, 1.249, 1.078, 0., 0.825, 0.)

HELIUM    = 2
NEON      = 10
ARGON     = 18
KRYPTN    = 36
XENON     = 54
RADON     = 86
########################


@auto_hash_equals struct MMFF94Parameters{T<:Real} <: AbstractForceFieldParameters
    sections::OrderedDict{String, BALLIniFileSection}
    radii::Vector{T}
    electronegativities::Vector{T}

    function MMFF94Parameters{T}(
            path::String = ball_data_path("forcefields/MMFF94/mmff94.ini")) where {T<:Real}

        ini_file = read_ball_ini_file(path, T)

        sections = ini_file.sections
        
        # see http://www.ccl.net/cca/data/MMFF94/
        radii = T.([
             0.33, 0.0,
             1.34, 0.90, 0.81, 0.77, 0.73, 0.72, 0.74, 0.0,
             1.54, 1.30, 1.22, 1.15, 1.09, 1.03, 1.01, 0.0,
             1.96, 1.74,
             1.44, 1.36, 0.00, 0.00, 0.00,
             0.00, 0.00, 0.00, 1.38, 1.31,
             1.19, 1.20, 1.20, 1.16, 1.15, 0.0,
             2.11, 1.92,
             1.62, 1.48, 0.00, 0.00, 0.00,
             0.00, 0.00, 0.00, 1.53, 1.48,
             1.46, 1.40, 1.41, 1.35, 1.33, 0.0
        ])
        
        electronegativities = T.([
             2.20, 0.0,
             0.97, 1.47, 2.01, 2.5, 3.07, 3.5, 4.10, 0.0,
             1.01, 1.23, 1.47, 1.74, 2.06, 2.44, 2.83, 0.0,
             0.91, 1.04,
             1.3, 1.5, 1.6, 1.6, 1.5,
             1.8, 1.8, 1.8, 1.9, 1.6,
             1.82, 2.02, 2.20, 2.48, 2.74, 0.0,
             0.89, 0.99,
             1.3, 1.4, 1.6, 1.8, 1.9,
             2.2, 2.2, 2.2, 1.9, 1.7,
             1.49, 1.72, 1.82, 2.01, 2.21, 0.0
        ])

        new(sections, radii, electronegativities)
    end
end

function MMFF94Parameters(path::String = ball_data_path("forcefields/MMFF94/mmff94.ini"))
    MMFF94Parameters{Float32}(path)
end


Base.show(io::IO, mmff_param::MMFF94Parameters{T}) where {T<:Real} = 
    print(io, 
        "MMFF94 parameters with $(length(mmff_param.sections)) sections.")


struct AtomTypeParams
    symbol::String7
    type::String3
    rule::String127
end
        
Base.show(io::IO, x::AtomTypeParams) = show(io::IO, (x.symbol, x.type, x.rule))

function _read_types_paramfile(path)
    lines = filter!(line -> !startswith(line, "*") ,readlines(path))
    lines = split.(lines, r"\s+\|\s+")
    types_param_dict = Dict{String3, Vector{AtomTypeParams}}()
    for line in lines
        vec = get!(types_param_dict, String3(line[1]), Vector{AtomTypeParams}())  
        push!(vec, AtomTypeParams(line[2:4]...))
    end
    return types_param_dict
end
        
        
        
function assign_atom_types_mmff94(ac::AbstractAtomContainer, mmff94_params, aro_rings)
    #add aro rings, add atomtype properties
    #read contents of TYPES.PAR
    path = ball_data_path("forcefields/MMFF94/TYPES.PAR")
    elem_to_params = _read_types_paramfile(path)

    ######################
    # setup aromatic
    ######################

    df = mmff94_params.sections["AtomTypeProperties"].data
    hetero_atom_types_df = @rsubset df begin
        :type != 10 && 
        (
            (:aspec == 7 && :crd == :val == 3) ||
            ((:aspec == 8 || :aspec == 16) && :crd == :val == 2)
        )
        @kwarg view = true
    end
    
    df = mmff94_params.sections["Aromatic"].data
    cation_df = @rsubset df begin
        !(endswith(:oldtype, "*"))
        :rsize == 5
        Bool(:imcat)
        @kwarg view = true
    end

    ring_5_df = @rsubset df begin
        !(endswith(:oldtype, "*"))
        :rsize == 5
        @kwarg view = true
    end
    
    
    map(mol -> assign_to(mol, elem_to_params), molecules(ac))

    ##############
    # assign heterocyclic 5 ring member types
    # e.g. CSA CSB N5A N5B C5
    ##############

    assign_heterocyclic_member_types(aro_rings,
        cation_df, 
        ring_5_df, 
        mmff94_params.sections["Symbols"].data,
        hetero_atom_types_df
    )


    ## assign hydrogens
    assign_hydrogen_types(ac,
        mmff94_params.sections["HydrogenTypes"].data,
        mmff94_params.sections["Symbols"].data
    )

    untyped_atoms_count = count(name -> name == "Any", atoms_df(ac).name)
    if untyped_atoms_count > 0
        @info "Could not assign atom type for $untyped_atoms_count atoms."
    end
end

function assign_hydrogen_types(ac, h_df, sym_df)
    hydro_adjacent_df = filter(at -> at.element == Elements.H || nbonds(at) == 1, atoms(ac))
    for at in hydro_adjacent_df
        nbonds(at) == 0 && continue
        partner = get_partner(first(bonds(at)), at)
        h_name = only(h_df[h_df.a .== partner.name, :b])
        H_nr = only(sym_df[sym_df.symbol .== h_name, :type])

        parse(Int8, at.atom_type) > H_nr && continue
        at.name = h_name
        at.atom_type = string(H_nr)
    end
end

function assign_heterocyclic_member_types(aro_rings, cation_df, ring_5_df, symbols_df, hetero_atom_types_df)
    for ring in filter(r -> length(r) == 6, aro_rings)
        for atom in ring
            if atom.name == "CNN+"
                atom.name = "CB"
                atom.atom_type = "37"
                map(bonds(atom)) do bond
                    p = get_partner(bond, atom)
                    p.atom_type == "55" && (p.atom_type = "40"; p.name = "NC=N")
                end
            end
        end
    end

    for ring in filter(r -> length(r) == 5, aro_rings)
        hetero_atom::Union{Missing, Atom} = missing
        cation_atom::Union{Missing, Atom} = missing
        anion_atom ::Union{Missing, Atom} = missing

        for at in ring
            if (
                    !in(at.element, (Elements.C, Elements.N)) || 
                    in(parse(Int, at.atom_type), hetero_atom_types_df.type)
                ) &&
                (
                    ismissing(hetero_atom) ||
                    hetero_atom.element == Elements.N
                )
                hetero_atom = at
            end

            at.name in cation_df.oldtype && (cation_atom = at)
            at.name in ("NM", "N5M") && (anion_atom = at)
        end

        charged = ismissing(hetero_atom) && any(ismissing, (cation_atom, anion_atom))
        charged_atom = ismissing(cation_atom) ? anion_atom   : cation_atom
        atom_x       = ismissing(hetero_atom) ? charged_atom : hetero_atom

        ismissing(atom_x) && continue
        L5 = 0

        for at in ring
            if hetero_atom == at        #test this ToDo
                L5 = 1
            elseif charged
                L5 = 4
            elseif is_bound_to(atom_x, at)
                L5 = 2
            else
                L5 = 3
            end
        end

        assign_5_ring_type(
            atom,
            L5,
            ismissing(anion_aom),
            ismissing(cation_atom),
            ring_5_df,
            symbols_df
        )

    end 
end
        
function assign_5_ring_type(atom, L5, is_anion, is_cation, ring_5_df, symbols_df)
    oldtype = atom.name
    is_anion && oldtype == "N5M" && return true
    #row = ring_5_df[(ring_5_df.oldtype .== oldtype) .& (ring_5_df.l5 .== L5), :]
    row = ring_5_df[findfirst(x-> x.oldtype == oldtype && x.l5 == L5, eachrow(ring_5_df))]

    found = true
    is_cation && !Bool(row.imcat)   && (found = false)
    is_anion  && !Bool(row.n5anion) && (found = false)

    if found 
        atom.name = row.aromtype
        atom.atom_type = string(symbols_df[findfirst(x -> x == row.arom_type, symbols_df.symbol), :type])
        return true
    end

    new_type = ""
    if atom.element == Elements.C
        L5 == 2 && (new_type = "C5A" )
        L5 == 3 && (new_type = "C5B" )
        L5 == 4 && (new_type = "C5"  )
    elseif atom.element == Elements.N
        L5 == 1 && (new_type = "NPYL")
        L5 == 2 && (new_type = "NSA" )
        L5 == 3 && (new_type = "N5B" )
        L5 == 4 && (new_type = "N5"  )

        oxide = false
        val   = 0

        for bond in bonds(atom)
            !in(bond.order, (BondOrder.unknown, BondOrder.Aromatic)) && (val += bond.order)
            bond.order == BondOrder.Aromatic && (val += 5)
            partner = get_partner(bond, atom)
            partner.element == Elements.O && nbonds(atom) == 1 && (oxide = true)
        end

        if val == 4
            if oxide
                new_type *= length(new_type) == 3 ? "X" : "OX"
            else
                new_type *= "+"
            end
        end

    elseif atom.element == Elements.O
        L5 == 1 && (new_type = "OFUR")
    elseif atom.element == Elements.S
        L5 == 1 && (new_type = "STHI")
    else
        @info "In \"mmff94_parameters.jl::assign_5_ring_type\" : No valid MMFF94 Type was assigned for $atom"
    end

    atom.name      = new_type
    atom.atom_type = string(symbols_df[findfirst(x -> x == new_type, symbols_df.symbol), :type])

end


function assign_to(mol, elem_to_params)
    #set some default type for all non-hydro atoms
    # untyped atoms will have atom_type = "0" and name = "Any"
    atoms_ = atoms_df(mol)
    atoms_.atom_type .= "0"         #default value for atom typer
    atoms_.name .= "Any"            #default name
    unassigned_atoms = Dict{String3, Set{Int}}()
    map(eachrow(atoms_)) do at
        push!(get!(unassigned_atoms, string(at.element), Set{Int}()), at.idx)
    end

    # each pair in elem_to_params amps an element to all SMARTS rules for that elem
    # assign the type of an atom if it is caught by one of the rules and
    # if it doesnt already have a type assigned
    for (elem, u_ats) in unassigned_atoms
        for tup in elem_to_params[elem]
            length(u_ats) == 0 && break
            
            m = @pipe SMARTSQuery(String(tup.rule)) |> match(_, mol) |> 
                atoms.(_) |>  Iterators.flatten(_)
            isempty(m) && continue
            untyped_atoms = (at for at in m if at.name == "Any" && string(at.element) == elem)
            
            # ignoring assignSpecificValues in atomTyper.C line 216
            map(untyped_atoms) do at
                at.name = tup.symbol
                at.atom_type = tup.type
                delete!(u_ats, at.idx)
            end
        end
    end

end


# `b` corresponds to the "BondStretch" file's column `b`, analogously for c
struct Stretch_Idx
    b::Int8
    c::Int8

    function Stretch_Idx(b::AbstractString, c::AbstractString)
        return b < c ? new(parse(Int8,b),parse(Int8,c)) : new(parse(Int8,c),parse(Int8,b))
    end

    function Stretch_Idx(b::Int, c::Int)
        return b < c ? new(b,c) : new(c,b)
    end
end

SI(b,c) = Stretch_Idx(b,c)

struct BondData{T<:Real}
    kb_normal::T
    r0_normal::T
    standard_bond_exists::Bool
                
    kb_sbmb::T
    r0_sbmb::T
    sbmb_exists::Bool
    empirical::Bool
end

#struct BondData?



function assign_bond_types_mmff94(ac::AbstractAtomContainer{T}, 
        mmff94_params,
        aromatic_rings) where T <: Real
    
    non_hydro_bonds = non_hydrogen_bonds(ac)

    ###### make stretch params
    stretch_params = Dict{Stretch_Idx, BondData{T}}()

    #may be super buggy
    map(eachrow(mmff94_params.sections["BondStretch"].data)) do row
        if Bool(row.a)
            stretch_params[SI(row.b, row.c)] = BondData(T(0), T(0), false, row.kb, row.r0, true, false)
        else
            stretch_params[SI(row.b, row.c)] = BondData(row.kb, row.r0, true, T(0), T(0), false, false)
        end
    end

    #######
    map(
        bond -> assign_bond_type(
                        bond,
                        stretch_params,
                        aromatic_rings,
                        mmff94_params.sections["AtomTypeProperties"].data
                    ),
        non_hydro_bonds
    )

    nothing
end


# assignMMFF94BondType, MMFF94.C l428
function assign_bond_type(bond, stretch_params, aromatic_rings, atom_types_df)

    # types 83-86 dont exist in the atoms type file.
    # this func allows indexing the atom type properties df
    at_idx(type) = parse(Int8, type) > 82 ? parse(Int8, type) - 4 : parse(Int8, type)

    # map a1.type + a2.type -> the entry in stretch_params
    a1_at = atom_by_idx(parent(bond), bond.a1)
    a2_at = atom_by_idx(parent(bond), bond.a2)
    bond_data = get(() -> missing, stretch_params, SI(a1_at.atom_type, a2_at.atom_type))
    ismissing(bond_data) && return false

    is_sbmb::Bool = bond_data.sbmb_exists && !bond_data.standard_bond_exists
    if bond_data.sbmb_exists == bond_data.standard_bond_exists
        if bond.order == BondOrder.Single

            #test this
            #potential problem: atom as from atoms() vs as from atoms_df
            at1_aromatic = a1_at in Iterators.flatten(aromatic_rings)
            at2_aromatic = a2_at in Iterators.flatten(aromatic_rings)

            if !(at1_aromatic || at2_aromatic) &&
                    Bool(atom_types_df[at_idx(a1_at.atom_type), :sbmb]) &&
                    Bool(atom_types_df[at_idx(a2_at.atom_type), :sbmb])
                is_sbmb = true
            elseif at1_aromatic && 
                    at2_aromatic &&
                    !any(ring -> a1_at in ring && a2_at in ring, aromatic_rings) 
                is_sbmb = true
            end
        end
    end
    is_sbmb && set_property!(bond, :MMFF94SBMB, true)
    !is_sbmb && haskey(bond.properties, :MMFF94SBMB) && (delete!(bond.properties, :MMFF94SBMB))

    return true

end


function get_partial_charge_(a,b,c, partial_charges_df, invalid_val, partial_bond_charges_df)
    b == c && return 0

    df = @rsubset partial_charges_df begin
        :a == a
        :b == b
        :c == c
        @kwarg view = true
    end
    
    isempty(df) && (@show empty;return invalid_val)
    only(df).charge != invalid_val && return only(df).charge

    #####heuristic value
    p1 = partial_bond_charges_df[partial_bond_charges_df.type .== b, :pbci]
    p2 = partial_bond_charges_df[partial_bond_charges_df.type .== c, :pbci]

    r = p2 - p1
    r == invalid_val && (@info "No ES parameters bt:$a at1:$b at2:$c.";)
    return r

end

#MMFF94Processors.C line 716
function assign_charges(ac::AbstractAtomContainer{T},
        mmff94_params, aro_rings) where T<:Real

    unassigned_atoms =  []
    invalid_val = 99.0
    # process "Charges"
        # rule_types : whereever :charge is "*". store :symbol
        # types_to_charges: :symbol => :charge
    charge_df = @rsubset mmff94_params.sections["Charges"].data begin
        :charge != "*"
        @kwarg view = true
    end

    rule_charge_df = @rsubset mmff94_params.sections["Charges"].data begin
        :charge == "*"
        @kwarg view = true
    end

    partial_charges_df = mmff94_params.sections["PartialCharges"].data

    # process "PartialBondCharges"
        # phis_ == :fcadj ?
    partial_bond_charges_df = mmff94_params.sections["PartialBondCharges"].data

    assign_partial_charges(ac, charge_df, rule_charge_df, aro_rings)
    charges = Dict{Atom{T}, T}()

    map(at -> setproperty!(at, :InitialCharge, at.charge), atoms(ac))
    for at in atoms(ac)
        nbonds(at) == 0 && continue
        phi = partial_bond_charges_df[parse(Int8, at.atom_type), :fcadj]
        phi == invalid_val && (push!(unassigned_atoms, at))
        phi == 0 && continue
        c = at.charge * phi
        at.charge = (1 - nbonds(at) * phi) * at.charge
        for bond in bonds(at)
            at2 = get_partner(bond, at)
            haskey(charges, at2) ? (charges[a2] += c) : (charges[a2] = c)

        end
    end

    #add all the charges up
    foreach((at, c) -> at.charge += c, charges)

    for at1 in atoms(ac)
        charge = copy(at1.charge)
        pcharge = 0

        for bond in bonds(at1)
            at2 = get_partner(bond, at1)
            at1_t = parse(Int, at1.atom_type)
            at2_t = parse(Int, at2.atom_type)
            bt = haskey(bond.properties, :MMFF94SBMB)

            if parse(Int8,at1.atom_type) < parse(Int8,at2.atom_type) 

                pcharge = (-1) * get_partial_charge_(
                    bt,at1_t,at2_t, partial_charges_df,
                    invalid_val, partial_bond_charges_df
                    )
            else
                pcharge = get_partial_charge_(
                    bt,at2_t,at1_t, partial_charges_df,
                    invalid_val, partial_bond_charges_df
                    )
            end
            abs(pcharge) == invalid_val && @show pcharge
            abs(pcharge) == invalid_val && (push!(unassigned_atoms, at1, at2); continue)
            charge += pcharge

        end

        at1.charge = charge
    end 

    return unassigned_atoms
end


function assign_partial_charges(
        ac::AbstractAtomContainer{T},
        charge_df,
        rule_charge_df,
        aro_rings) where T<:Real


    atoms_df(ac).charge .= 0
    # get all atoms for which the charge can just be read from charge_df
    atoms_simple_charge = @rsubset atoms_df(ac) begin
        !(in(:name, rule_charge_df.symbol))
        @kwarg view = true
    end

    for at in eachrow(atoms_simple_charge)
        idx = findfirst(x -> x == at.name, charge_df.symbol)
        isnothing(idx) && continue
        at.charge = parse(T,charge_df[idx, :charge])
    end

    #these can push atoms into N5Ms and NIMs
    #do N5M find all arom rings which have atom
    
    aro_rings_ = [isempty(ring) ? [] : first(ring) for ring in aro_rings]
    N5Ms = [[ at for at in ring if at.name == "N5M" ] for ring in aro_rings_]
    NIMs = [[ at for at in ring if at.name == "Nim+"] for ring in aro_rings_]

    #do O2S
    atoms_O2S = @rsubset atoms_df(ac) begin
        :name == "O2S"
        @kwarg view = true
    end

    for atom in eachrow(atoms_O2S)
        for bond in bonds(atom)
            partner = get_partner(bond, atom)
            if partner.element == Elements.S && nbonds(partner) == 3
                nr = 0
                map(bonds(partner)) do p_bond
                    !in(p_bond.order, (BondOrder.unknown, BondOrder.Aromatic)) && 
                        (nr += p_bond.order)
                    p_bond.order == BondOrder.Aromatic && 
                        (nr += 5)
                end
                nr == 3 && (atom.charge = -0.5; break)
            end
        end
    end

    #do NR%
    atoms_NR = @rsubset atoms_df(ac) begin
        :name == "NR%"
        @kwarg view = true
    end

    for atom in eachrow(atoms_NR)
        for bond in bonds(atom)
            bond.order == BondOrder.Triple                    && 
                get_partner(bond, atom).element == Elements.N &&
                (atom.charge = 1;)
        end
    end

    #do SM
    atoms_NR = @rsubset atoms_df(ac) begin
        :name == "SM"
        @kwarg view = true
    end

    for atom in eachrow(atoms_NR)
        partner = get_partner(only(bonds(atom)), atom)
        for bond in bonds(partner)
            partner2 = get_partner(bond, partner)
            if partner2.element == Elements.O && nbonds(partner2) == 1 &&
                only(bonds(partner2)).order == BondOrder.Single
                atom.charge = -0,5
            end
        end

        arom.charge = -1
    end

    #handle nim, n5m
    for ring in (r for r in N5Ms if !isempty(r))
        charge_ = -1/length(ring)
        map(at -> at.charge = charge_, ring)
    end

    for ring in (r for r in NIMs if !isempty(r))
        nr_partner = 0
        for at in ring
            for bond in bonds(ring)
                get_partner(bond, at).name in ("NGD+", "NCN+") && (nr_partner +=1) 
            end
        end

        charge_ = 1.0 / (length(ring) + nr_partner)
        map(at -> at.charge = charge_, ring)
    end

end