@testitem "aromaticity" begin 

    result = [0, 6, 0, 6, 6, 6, 0, 6, 0]

    system = load_sdfile(ball_data_path("../test/data/descriptors_test.sdf"))
    rings = map(find_sssr, molecules(system))
    aro_rings = map(r -> map(aromatize_simple,r), rings)

    bond_counts = map(mol -> count(bond -> bond.order == BondOrder.Aromatic, bonds(mol)), molecules(system))
    @test all(result .== bond_counts)

    at_counts = map(mol -> count(at -> has_property(at, :IsAromatic), atoms(mol)), molecules(system))
    @test all(result .== at_counts)


end