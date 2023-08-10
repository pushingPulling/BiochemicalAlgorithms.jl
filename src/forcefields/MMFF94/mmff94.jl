export MMFF94FF, _read_types_paramfile

function _check_atom_type(atom; unassigned_atoms)
    is_empty(atom.atom_type) && atom in unassigned_atoms && return false
    push!(unassigned_atoms, atom)
    @info "invalid atom type"
    return false
end

get_mmff94_default_options(T=Float32) = Dict{Symbol, Any}(
    # copied from BALL
    :nonbonded_cutoff               => T(16.0),
    :vdw_cutoff                     => T(15.0),
    :vdw_cuton                      => T(13.0),
	:electrostatic_cutoff           => T(15.0),
    :electrostatic_cuton            => T(13.0),
    :distance_dependent_dielectric  => false,
    :assign_charges                 => true,
    :assign_typenames               => true,
    :assign_types                   => true,
    :overwrite_nonzero_charges      => true,
    :overwrite_typenames            => true,
    :periodic_boundary_conditions   => false,       #missing in BAll
    :periodic_box_width             => T(100.0),    #missing in BALL
    :periodic_box_height            => T(100.0),    #missing in BALL
    :periodic_box_depth             => T(100.0),    #missing in BALL
    :max_number_of_unassigned_atoms => typemax(Int32),   #missing in BALL
    :MMFF_ES_ENABLED                => true,
    :MMFF_VDW_ENABLED               => true
)

function MMFF94FF(
        ac::AbstractAtomContainer{T}, 
        filename=ball_data_path("forcefields/MMFF94/mmff94.ini");   # ToDo: What does this file have??
        constrained_atoms=Vector{Int}(),
        options = get_mmff94_default_options(T)) where {T<:Real}            # ToDo: do i want this?

    mmff94_params = MMFF94Parameters(filename)                      
    mmff94_ff = ForceField{T}(                                      # ToDo: implement this function call
        "MMFF94",
        ac, 
        mmff94_params, 
        options,
        Dict{String, AtomTypeTemplate{T}}(),                          #ToDO: Find this function and verify its usage
        Vector{AbstractForceFieldComponent{T}}(),
        Dict{String, T}(),
        Vector{Atom{T}}(),
        constrained_atoms
    )


    all_rings = map(find_sssr, molecules(ac))
    #for some reason map splats the returned list if length(molecules(ac)) == 1 ...
    #length(molecules(ac)) == 1 && (all_rings = only(all_rings);)
    
    #map(mol_rings -> map(ring -> map(at -> set_property!(at, :in_ring, true), ring), mol_rings), all_rings)
    #for mols_of_sys in all_rings
    map(rings_of_mol -> map(ring -> map(atom -> set_property!(atom, :in_ring, true), ring), rings_of_mol), all_rings)
        #for rings_of_mol in all_rings
        #    for ring in rings_of_mol
        #        for atom in ring
        #            set_property!(atom, :in_ring, true)
        #        end
        #    end
        #end
    #end

    aromatic_rings = map(aromatize_simple, all_rings)


    unassigned_bonds = map(kekulizer!, molecules(ac), aromatic_rings) 

    length(collect(Iterators.flatten(unassigned_bonds))) > 0 && @info "There are $(length(collect(Iterators.flatten(unassigned_bonds)))) unassigned Bonds in mmff94.jl::MMFF94FF."
    #fix? #length(unassigned_atoms) > 0 && push!(mmff94_ff.unassigned_atoms, Iterators.flatten(unassigned_atoms)...)

    if mmff94_ff.options[:assign_types]    
        assign_atom_types_mmff94(ac, mmff94_params, aromatic_rings)
    end

    assign_bond_types_mmff94(ac, mmff94_params, all_rings)

    if mmff94_ff.options[:assign_charges]
        unassigned_atoms = assign_charges(ac, mmff94_params, aromatic_rings)
        length(unassigned_atoms) > 0 && @info "Could not assign charges for $(length(unassigned_atoms)) atoms."
        length(unassigned_atoms) > 0 && push!(mmff94_ff.unassigned_atoms, unassigned_atoms...)
    end


    append!(            # ToDo: implement all of these
        mmff94_ff.components,
        [
            MStretchBendComponent{T}(mmff94_ff, all_rings, aromatic_rings),
            MTorsionComponent{T}(mmff94_ff, all_rings, aromatic_rings),
            MNonBondedComponent{T}(mmff94_ff),
            MOutOfPlaneComponent{T}(mmff94_ff)
        ]
    )

    mmff94_ff
end


