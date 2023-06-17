export MMFF94FF

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
    :max_number_of_unassigned_atoms => typemax(Int32)   #missing in BALL
)

function MMFF94FF(
        ac::AbstractAtomContainer{T}, 
        filename=ball_data_path("forcefields/MMFF94/mmff94.ini");   # ToDo: What does this file have??
        constrained_atoms=Vector{Int}()) where {T<:Real}            # ToDo: do i want this?

    mmff94_params = MMFF94Parameters(filename)                      
    mmff94_ff = ForceField{T}(                                      # ToDo: implement this function call
        "MMFF94",
        ac, 
        mmff94_params, 
        get_mmff94_default_options(T),
        init_atom_types(mmff94_params, T),                          #ToDO: Find this function and verify its usage
        Vector{AbstractForceFieldComponent{T}}(),
        Dict{String, T}(),
        Vector{Atom{T}}(),
        constrained_atoms
    )


    ######Overall program#######
    #   insert componentes
    #   specific setup
    ###########################

    ##### specific setup MMFF94.C@line146

    #bonds that arent hydrogen bonds
    #(BALL uses the concept of aromatized rings - which rings get aromatized?)
    non_hydro_bonds = non_hydrogen_bonds(ac)
    
    all_rings = map(find_sssr, molecules(ac))

    ###-###
    # Todo: make function that aromatizes a ring
    aromatic_rings = map(aromatize_simple, all_rings)

    ###-###



    non_aromatic_rings = setdiff(all_rings, aromatic_rings)

    #simply collect all bonds which are in an aromatic ring
    aromatic_bonds = filter(b -> atom_by_idx(parent_system(mol), b.a1) in ring && atom_by_idx(parent_system(mol), b.a2) in ring, bonds(mol))
    # get all unassigned bonds, throw error for each bond that cant be kekulized
    unassigned_bonds = nothing 

    # => kekuliser_.setup for each molecule

    # if assign types - lots of stuff happens here
    # unconditional assign bondtypes - lots happens here
    # if assign_charges - 

    ##### specific setup end


    append!(                                                        # ToDo: implement all of these
        mmff94_ff.components,
        [
            #QuadraticStretchComponent{T}(mmff94_ff),
            #QuadraticBendComponent{T}(mmff94_ff),
            #TorsionComponent{T}(mmff94_ff),
            #NonBondedComponent{T}(mmff94_ff)
        ]
    )

    mmff94_ff
end