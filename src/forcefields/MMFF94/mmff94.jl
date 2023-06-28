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
        Dict{String, AtomTypeTemplate{T}}(),                          #ToDO: Find this function and verify its usage
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

    #bonds that arent hydrogen bonds)
    
    all_rings = map(find_sssr, molecules(ac))
    aromatic_rings = map(aromatize_simple, all_rings)

    #TODO: may be a bit buggy, need to test
    # see if this needs to be flattened
    non_aromatic_rings = setdiff(all_rings, aromatic_rings)

    #simply collect all bonds which are in an aromatic ring
    #assume that bond order Aromatic is set, due to aromatize_simple
    aromatic_bonds = filter(b -> b.order == BondOrder.Aromatic,bonds_df(ac))
    # => kekuliser_.setup for each molecule
    # get all unassigned bonds, throw error for each bond that cant be kekulized
    #TODO: see if that list needs to be flattened
    unassigned_atoms = map(kekulizer!, molecules(ac), aromatic_rings) 
    @info "There are $(length(collect(Iterators.flatten(unassigned_atoms)))) unassigned Bonds in mmff94.jl::MMFF94FF."

    #initialise params ll 230-228
    #already done in mmff94_params
    #ll. 243-248 see below
    #ll. 249,250 see below

    ######################################################
    # set atom types
    ######################################################

    if mmff94_ff.options[:assign_types]
        assign_atom_types_mmff94(ac, mmff94_params, aromatic_rings)
    end

    ######################################################
    # set bond types
    ######################################################
    assign_bond_types_mmff94(ac, mmff94_params, all_rings)
    ######################################################
    # ToDo: set charges
    ######################################################

    if mmff94_ff.options[:assign_charges]
        #setup charges processor
        #set ESparameters
        #types_to_charges_[key]:: str -> float ... symbol -> charge
        #rule_types set of str. come from setup where "*" is a key
        #read partial charges and partial bond charges

        #apply charge processor l716
        unassigned_atoms = assign_charges(ac, mmff94_params, aromatic_rings)
        length(unassigned_atoms) > 0 && @info "Could not assign charges for $(length(ua)) atoms."
            
        #set unassigned atoms for later processing
    end


    append!(            # ToDo: implement all of these
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
