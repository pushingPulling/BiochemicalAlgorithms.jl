@testitem "kekulizer" begin 

    system = load_sdfile(ball_data_path("../test/data/descriptors_test.sdf"))
    rings = map(find_sssr, molecules(system))

    c = map(mol -> count(bond -> bond.order == BondOrder.Aromatic, bonds(mol)), molecules(system))
    @test c == zeros(9)

    aro_rings = aromatize_simple.(rings)

    c = map(mol -> count(bond -> bond.order == BondOrder.Aromatic, bonds(mol)), molecules(system))
    @test c == [0, 6, 0, 6, 6, 6, 0, 6, 0]

    kekulizer!.(molecules(system), aro_rings)

    c = map(mol -> count(bond -> bond.order == BondOrder.Aromatic, bonds(mol)), molecules(system))
    @test c == zeros(9)
end