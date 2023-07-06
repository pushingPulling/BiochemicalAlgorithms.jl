using BiochemicalAlgorithms
#using MolecularGraph

function check_sssr()
    system = load_sdfile(ball_data_path("../test/data/rings_test.sdf"))

    num_ring_atoms = [6, 6, 6, 6, 6, 13, 0, 6, 4]

    detected_num_ring_atoms = map(sum, map(is_ring_atom, molecules(system)))

    num_sssrs      = [1, 1, 1, 1, 1, 3, 0, 1, 3]
    num_sssr_atoms = [[6], [6], [6], [6], [6], [6, 3, 6], [], [6], [3, 3, 3]]

    detected_sssrs = map(find_sssr, molecules(system))
end

function check_aro()

    system = load_sdfile(ball_data_path("../test/data/descriptors_test.sdf"))
    mol = molecules(system)[2]

    ring = find_sssr(mol)[1]
    bonds_in_ring = Set()
    map(l -> push!(bonds_in_ring, bonds(l)...), ring)
    orders = [bond.order for bond in bonds_in_ring]

    a = [deepcopy(at) for at in ring]
    aro_rings = aromatize_simple(ring)
    b = [at for at in ring] 
end

function test_smart()
    sys = System()
    mol = Molecule(sys)
    BiochemicalAlgorithms.Atom(mol, 2, Elements.O)
    BiochemicalAlgorithms.Atom(mol, 3, Elements.C)
    BiochemicalAlgorithms.Atom(mol, 4, Elements.O)
    BiochemicalAlgorithms.Atom(mol, 5, Elements.H)
    #BiochemicalAlgorithms.Atom(mol, 6, Elements.H)

    BiochemicalAlgorithms.Bond(sys, 2,3, BondOrder.Double)
    BiochemicalAlgorithms.Bond(sys, 3,4, BondOrder.Single)
    BiochemicalAlgorithms.Bond(sys, 4,5, BondOrder.Single)
    #BiochemicalAlgorithms.Bond(sys, 5,6, BondOrder.Single)
    println(bonds(mol)[1])

    c = SMARTSQuery("[CX3](=O)[OX2H1]")

    m = match(c, mol)
end

function test_kek()
    system = load_sdfile(ball_data_path("../test/data/descriptors_test.sdf"))
    rings = map(find_sssr, molecules(system))


    c = map(mol -> count(bond -> bond.order == BondOrder.Aromatic, bonds(mol)), molecules(system))
    @show "pre arom", c

    aro_rings = aromatize_simple.(rings)

    c = map(mol -> count(bond -> bond.order == BondOrder.Aromatic, bonds(mol)), molecules(system))
    @show "post arom", c

    kekulizer!.(molecules(system), aro_rings)

    c = map(mol -> count(bond -> bond.order == BondOrder.Aromatic, bonds(mol)), molecules(system))
    @show "post kek", c
end

function test_read_types()
    system = load_sdfile(ball_data_path("../test/data/descriptors_test.sdf"))
    mol = molecules(system)[2]
    @show count( at -> at.element == Elements.C, atoms(mol))
    @show count( at -> at.element == Elements.C && at.atom_type == "-1", atoms(mol))

    path = ball_data_path("forcefields/MMFF94/TYPES.PAR")
    elem_to_params = _read_types_paramfile(path)

    assign_to(elem_to_params, mol)

    @show [at.atom_type for at in atoms(mol)]
    @show count( at -> at.atom_type == "-1", atoms(mol))
end


function test_mm()
    i = 1
    system = load_sdfile(ball_data_path("../test/data/descriptors_test.sdf"))
    mmff94 = MMFF94FF(system)

    #using Serialization
    #serialize("zff", mmff94)
    #serialize("zsys", system)
    #mmff94 = deserialize("zff")
    #sytem = deserialize("zsys")

    msb = BiochemicalAlgorithms.MStretchBendComponent{Float32}(mmff94, "msbc", Dict{Symbol, Any}(),
        Dict{String, Float32}(), Dict{NTuple{2,Atom{Float32}},BiochemicalAlgorithms.MStretch}(), 
        BiochemicalAlgorithms.MBend[], BiochemicalAlgorithms.MStretchBend[])
    #msb = deserialize("zmsb")
    #serialize("zmsb", msb)
    #all_rings = deserialize("zar")
    #aromatic_rings = deserialize("zaro" )
    all_rings = map(find_sssr, molecules(system))
    #serialize("zar", all_rings)
    aromatic_rings = map(aromatize_simple, all_rings)
    #serialize("zaro", aromatic_rings)
    setup!(msb, aromatic_rings, all_rings)
end
using Profile
using StatProfilerHTML
test_mm()
Profile.clear_malloc_data()
@time test_mm()
#@profilehtml test_mm()


#(c, mol.idx) = (12, 1)
#(c, mol.idx) = (6, 38)
#(c, mol.idx) = (15, 63)
#(c, mol.idx) = (12, 112)
#(c, mol.idx) = (13, 155)
#(c, mol.idx) = (16, 200)
#(c, mol.idx) = (1, 263)
#(c, mol.idx) = (6, 273)
#(c, mol.idx) = (0, 298)


# build ------>stretchbend<----------------
    # consists of a bend and a stretch
    # stretchbend is an emergent quality of those two

# find out how i make this an AbstractForceFieldComponent
    # see if there isa  template which i should use
        # yes -> copy funcitons and fill in with specific type
    # learn about mmff94 stretchbend more: Which params do i need etc?

# read params
    # this time, figure out how ball accesses thigns and how i should access things
        # AngleBend
        # StretchBend
        # StretchBendEmpirical
        # equivalencies?
        # BondStretch?          -> from bond_parameters in MMFF94.C
        # EmpiricalBondParams?  -> from bond_parameters in MMFF94.C

# lots of setup
    # learn the theory behind the values in detail
    # what they stand for
    # where the values are from
    # intended behaviour; physical mechanics

# learn unitful

# learn the actual mmff94 formulas with units
    # N = kg * m / s^2
    # J = N * m
    # kcal = J    
    # md = Joule/A = N

    # is morse function 4th order, returns kcal/mol

# make my own stretchbend class
