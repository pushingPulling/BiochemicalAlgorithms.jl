using BiochemicalAlgorithms
using MolecularGraph

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

