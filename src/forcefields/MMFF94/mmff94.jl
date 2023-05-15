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

    # ToDo: Do all needed setups here. Might be a lot of work to figure things out.
    # maybe check what c++ BALL does and how the data looks there. Compare to how data looks here

    # find all rings using sssr?
    #solution : find_sssr()


    #Todo: Apparently kekulize all rings and all aromatic rings
    #todo: What exactly is kekulize lmao - finding alternative conformation of molecules by shuffling the double bonds around

    #Todo: decide on a solution
    # possibly Moleculargraph kekulizer? But does it even work?
    # if yes: find a way to convert a MolGraph back to a molecule
    # if no: copy ball's kekulizer and use Molecular Graphs smart matcher


    #Todo: Read param ll 230-250 in MMFF94.C

    #todo: assign types
    #todo assign bond types




    assign_typenames_and_charges!(mmff94_ff)                         # ToDo: see if i want this

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