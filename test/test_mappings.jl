using Statistics: mean
using LinearAlgebra: Hermitian, eigen
@testset "Atom_Bijection:" begin
    mol = mol = load_pubchem_json("data/aspirin_pug.json")
    mol2 = deepcopy(mol) 

    ab = TrivialAtomBijection(mol, mol2)
    @test ab isa TrivialAtomBijection
    @test size(ab.atoms_A) == size(ab.atoms_B)

end

@testset "Rigid_Mapping: RigidTransform" begin
    v = Vector3{Float32}(1, 1, 1)
    m = Matrix3{Float32}(1, 2, 3, 4, 5, 6, 7, 8, 9)
    r = RigidTransform(m,v)

    @test r isa RigidTransform
    @test r.rotation isa Matrix3{Float32}
    @test r.translation isa Vector3{Float32}
end

@testset "Rigid_Mapping: Function translate" begin
    
    mol = mol = load_pubchem_json("data/aspirin_pug.json")
    mol2 = deepcopy(mol) 

    v = Vector3{Float32}(1.0,1.0,1.0)
    translate!(mol2, v)

    for i in eachindex(mol2.atoms.r)
        @test isapprox(mol2.atoms.r[i][1], mol.atoms.r[i][1] + 1.0)
        @test isapprox(mol2.atoms.r[i][2], mol.atoms.r[i][2] + 1.0)
        @test isapprox(mol2.atoms.r[i][3], mol.atoms.r[i][3] + 1.0)
    end
end

@testset "Rigid_Mapping: Function rigid_transform" begin
  
    mol = Molecule()
    # add atoms
    for i in 1:3
        atom = ( number = i,
        name = "my fancy atom", 
        element = Elements.H, 
        atomtype = "heavy",
        r = Vector3{Float64}(i*1.0, i*2.0, 0.0),
        v = Vector3{Float64}(0.0, 0.0, 0.0),
        F = Vector3{Float64}(0.0, 0.0, 0.0),
        has_velocity = false,
        has_force = false,
        frame_id = 1,
        properties = Dict{String, Any}()
        )
        push!(mol, atom)
        @test count_atoms(mol) == i
    end

    mol2 = deepcopy(mol) 
    
    # performs counter clockwise rotation by 90 degree with translation by v=(0, 0, 0)
    v = Vector3{Float32}(0, 0, 0)
    m = Matrix3{Float32}(1, 0, 0, 0, 0, -1, 0, 1, 0)
    r = RigidTransform(m,v)

    rigid_transform!(mol2, r)

    # check first atom
    @test isapprox(mol2.atoms.r[1][1], 1.0)
    @test isapprox(mol2.atoms.r[1][2], 0.0)
    @test isapprox(mol2.atoms.r[1][3], -2.0)

    # check second atom
    @test isapprox(mol2.atoms.r[2][1], 2.0)
    @test isapprox(mol2.atoms.r[2][2], 0.0)
    @test isapprox(mol2.atoms.r[2][3], -4.0)

    # check third atom
    @test isapprox(mol2.atoms.r[3][1], 3.0)
    @test isapprox(mol2.atoms.r[3][2], 0.0)
    @test isapprox(mol2.atoms.r[3][3], -6.0)

    # performs counter clockwise rotation by 90 degree with translation by v=(1, 1, 1)
    mol2 = deepcopy(mol) 
    v = Vector3{Float32}(1, 1, 1)
    m = Matrix3{Float32}(1, 0, 0, 0, 0, -1, 0, 1, 0)
    r = RigidTransform(m,v)

    rigid_transform!(mol2, r)

    # check first atom
    @test isapprox(mol2.atoms.r[1][1], 2.0)
    @test isapprox(mol2.atoms.r[1][2], 1.0)
    @test isapprox(mol2.atoms.r[1][3], -1.0)

    # check first atom
    @test isapprox(mol2.atoms.r[2][1], 3.0)
    @test isapprox(mol2.atoms.r[2][2], 1.0)
    @test isapprox(mol2.atoms.r[2][3], -3.0)

    # check first atom
    @test isapprox(mol2.atoms.r[3][1], 4.0)
    @test isapprox(mol2.atoms.r[3][2], 1.0)
    @test isapprox(mol2.atoms.r[3][3], -5.0)
    
end

@testset "Rigid_Mapping: Function compute_rmsd" begin
    mol = mol = load_pubchem_json("data/aspirin_pug.json")
    
    ab = TrivialAtomBijection(mol, mol)
    @test isapprox(compute_rmsd(ab), 0.0)

    mol2 = deepcopy(mol) 
    translate!(mol2, Vector3{Float32}(1.0, 0, 0))
    ab2 = TrivialAtomBijection(mol,mol2)
    @test isapprox(compute_rmsd(ab2), 1.0)

end

@testset "Rigid_Mapping: Function compute_rmsd_minimizer" begin
    mol = Molecule{Float64}()
    # add atoms
    for i in 1:3
        atom = ( number = i,
        name = "my fancy atom", 
        element = Elements.H, 
        atomtype = "heavy",
        r = Vector3{Float64}(i*1.0, i*2.0, 0.0),
        v = Vector3{Float64}(0.0, 0.0, 0.0),
        F = Vector3{Float64}(0.0, 0.0, 0.0),
        has_velocity = false,
        has_force = false,
        frame_id = 1,
        properties = Dict{String, Any}()
        )
        push!(mol, atom)
        @test count_atoms(mol) == i
    end

    # first: consider simple translation of coordinates by a vector v=(1,1,0)
    mol2 = deepcopy(mol) 
    mol2.atoms.r[1] = Vector3{Float64}(2.0, 3.0, 0.0)
    mol2.atoms.r[2] = Vector3{Float64}(3.0, 5.0, 0.0)
    mol2.atoms.r[3] = Vector3{Float64}(4.0, 7.0, 0.0)

    ab = TrivialAtomBijection{Float64}(mol, mol2)
    rt = compute_rmsd_minimizer(ab)
    @test rt isa RigidTransform
    @test rt.rotation isa Matrix3{Float64}
    @test rt.translation isa Vector3{Float64}
    @test rt.translation == Vector3{Float64}(1.0, 1.0, 0.0)

    # second: consider simple rotation of coordinates without translation
    mol2.atoms.r[1] = Vector3{Float64}(1.0, 0.0, -2.0)
    mol2.atoms.r[2] = Vector3{Float64}(2.0, 0.0, -4.0)
    mol2.atoms.r[3] = Vector3{Float64}(3.0, 0.0, -6.0)

 

    r_A = mol.atoms.r
    r_B = mol2.atoms.r
    mean_A = mean(r_A)
    mean_B = mean(r_B)


    R = mapreduce(t -> t[1] * transpose(t[2]), +, zip(r_B .- Ref(mean_B), r_A .- Ref(mean_A)))

    C = Hermitian(transpose(R) * R)
  
    μ, a = eigen(C)

    println("mu ",μ )
    println("a", a)

    ab2 = TrivialAtomBijection{Float64}(mol, mol2)
    rt = compute_rmsd_minimizer(ab2)
    println(rt.translation)
    println(rt.rotation)
end
