using BiochemicalAlgorithms
using Unitful
using LinearAlgebra
BA = BiochemicalAlgorithms
#using MolecularGraph
using Serialization
using BenchmarkTools


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


function test_mm_ser()
    syst = load_sdfile(ball_data_path("../test/data/descriptors_test_m8.sdf"))
    mmff94 = MMFF94FF(syst)

    serialize("zff", mmff94)
    serialize("zsys", syst)
end

function test_mm()
    system = load_sdfile(ball_data_path("../test/data/descriptors_test_m8.sdf"))
    mmff94 = MMFF94FF(system)
    #compute_energy(mmff94)
    #compute_forces(mmff94)
end

function test_mm_deser()
    mmff94 = deserialize("zff")
    system = deserialize("zsys")
    ring   = deserialize("zring")
    all_rings   = deserialize("zrings")

    return mmff94, system, ring, all_rings
end

function test()
    #using LinearAlgebra
    #using Unitful
    #const BA = BiochemicalAlgorithms
    charmm_force_factor = (1u"cal" |> u"J" |> ustrip) * (1000 * 1e10 / ustrip(Unitful.Na))
    #
    syst = load_sdfile(ball_data_path("../test/data/MMFF94-bend.sdf"))
    mm = MMFF94FF(syst)
    #@test natoms(syst) == 3
    a1, a2, a3 = atoms(syst)
    @show charmm_force_factor
    @show a1.atom_type, a2.atom_type, a3.atom_type
    #
    mm.options[:bends_enabled] = true
    mm.options[:stretches_enabled] = false
    mm.options[:stretch_bends_enabled] = false
    #
    sb = mm.components[1]
    @show length(sb.bends), length(sb.stretch_bends), length(sb.stretches)
    bend = sb.bends[1]
    compute_forces(sb)
    energy = compute_energy(sb)
    @show "nnnnsss", bend.n1, bend.n2
    @show bend.theta0, bend.theta_d, bend.theta, bend.ka, bend.is_linear, bend.ATIJK
    precision = 0.5
    charmm_energy = 3.09341u"cal" |> u"J" |> ustrip
    @show energy, charmm_energy
    @show abs(energy - charmm_energy), precision
    v1 = Vector3(0          , 27.3344889 , 0) * -(charmm_force_factor)
    v2 = Vector3(27.3344889 , -27.3344889, 0) * -(charmm_force_factor)
    v3 = Vector3(-27.3344889, 0          , 0) * -(charmm_force_factor)
    #
    precision = 2e-12
    @show a1.F, v1, a1.F .== v1
    @show a2.F, v2, a2.F .== v2
    @show a3.F, v3, a3.F .== v3
    #
    @show abs(v2.y - a2.F.y), precision
    #
    a3.r = (0., 2.96900, 0)
    v1 = Vector3(0          , 27.3344889 , 0) * -(charmm_force_factor)
    v2 = Vector3(8.92122591 , -27.3344889, 0) * -(charmm_force_factor)
    v3 = Vector3(-8.92122591, 0          , 0) * -(charmm_force_factor)
    #
    compute_forces(sb)
    #
    @show a1.F, v1, a1.F .== v1
    @show a2.F, v2, a2.F .== v2
    @show a3.F, v3, a3.F .== v3
    #
    @show abs(v2.y - a2.F.y)
    #
end

using Profile
using StatProfilerHTML
using DataFrames
using DataFramesMeta
#test_mm_ser()
#mm, syst = test_mm_deser()
#tors = mm.components[2]
#moop = mm.components[4]
#setup!(moop)
#println(length(moop.out_of_plane_bends))

using LinearAlgebra
using Unitful
using Test
joule_per_cal = 4.184
charmm_force_factor = (1u"cal" |> u"J" |> ustrip) * (1000 * 1e10 / ustrip(Unitful.Na))
force_prefactor = ustrip(u"kJ/mol/angstrom"/Unitful.Na |> u"N")

syst = load_sdfile(ball_data_path("../test/data/MMFF94-torsion.sdf"))
mm = MMFF94FF(syst)
@test natoms(syst) == 4
a1, a2, a3, a4 = atoms(syst)

tc = mm.components[2]
@test length(tc.torsions) == 1

v1 = Vector3( 0             ,   0,  12) * -(charmm_force_factor)
v2 = Vector3(-5.19556474E-16,   0, -12) * -(charmm_force_factor)
v3 = Vector3(-6             ,   0, -6)  * -(charmm_force_factor)
v4 = Vector3( 6             ,   0,  6)  * -(charmm_force_factor)

compute_forces(tc)
precision = 1e-16

a1, a2, a3, a4 = atoms(syst)
@test all(isapprox.(a1.F, v1, atol=precision))
@test all(isapprox.(a2.F, v2, atol=precision))
@test all(isapprox.(a3.F, v3, atol=precision))
@test all(isapprox.(a4.F, v4, atol=precision))

charmm_energy = 6.0u"cal" |> u"J" |> ustrip
@test isapprox(compute_energy(tc), charmm_energy, atol=0.1)
@test abs(compute_energy(tc) - charmm_energy) < 1e-5


#=
root = ball_data_path("../test/data/")
names = ["descriptors_test.sdf", "MMFF94_test1.sdf", "MMFF94_test2.sdf", "MMFF94-bend-lin.sdf", "MMFF94-bend.sdf", "MMFF94-bend2.sdf", "MMFF94-bend3.sdf", "MMFF94-plane.sdf", "MMFF94-plane2.sdf", "MMFF94-stretch.sdf", "MMFF94-torsion.sdf", "MMFF94-vdw.sdf", "MMFF94-vdw2.sdf", "rings_test.sdf", "sdfile_test_1.sdf"]
for name in names
    syst = load_sdfile(root*name)
    mm = MMFF94FF(syst)
    sbc = mm.components[1]
    tors = mm.components[2]
    nbc = mm.components[3]
    moop = mm.components[4]

    println("name $(name)")
    println("moops: $(length(moop.out_of_plane_bends));  tors: $(length(tors.torsions)); interacts: $(length(nbc.es_interactions)); stretchBends: $(length(sbc.stretch_bends)); streches $(length(sbc.stretches)); bends $(length(sbc.bends))\n")
end
=#

#just after assign type and before assign heterocyc type
#=
atoms_df(syst).atom_type .= ["2","2","54","57","62","2","39","2","9","40","40","0","22","22","22","2","2","1","1","1","1","6","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0"]
atoms_df(syst).name .= ["CSP2","CSP2","N+=C","CGD+","NM","CSP2","NPYL","CSP2","N=C","NC=C","NC=N","Any","CR3R","CR3R","CR3R","C=C","C=C","CR","CR","CR","CR","OR","Any","Any","Any","Any","Any","Any","Any","Any","Any","Any","Any","Any","Any","Any","Any","Any","Any"]

aro_rings = mm.components[1].rings_cache[2]

df = mm.parameters.sections["AtomTypeProperties"].data
hetero_atom_types_df = @rsubset df begin
    :type != 10 && 
    (
        (:aspec == 7 && :crd == :val == 3) ||
        ((:aspec == 8 || :aspec == 16) && :crd == :val == 2)
    )
    @kwarg view = true
end

df = mm.parameters.sections["Aromatic"].data
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


BA.assign_heterocyclic_member_types(aro_rings,
    cation_df, 
    ring_5_df, 
    mm.parameters.sections["Symbols"].data,
    hetero_atom_types_df
)
=#



#map(aromatize_simple, all_rings)
#aromatize_simple(ring)




#a1, a2 = eachrow(atoms_df(syst))

    #C++ BALL has an error where it assigns a negative charge to brom when it shouldnt
    # it also assigns Zinc a +2 chage when it shouldnt
    # reimplementing this error so its consistent to old Test
    #mm = MMFF94FF(syst)

 #   a_brom = only(filter(a -> a.element == Elements.Br, (a1,a2)))
 #   a_brom.charge = -1
 #   a_brom.atom_type = "91"
 #   a_brom.name = "BR-"
 #   a2.charge = 2
 #   a2.atom_type = "95"
 #   a2.name = "ZN+2"
#
 #   mm.options[:electrostatic_cuton] = 8
 #   mm.options[:electrostatic_cutoff] = 12
 #   mm.options[:MMFF_ES_ENABLED] = true
 #   mm.options[:MMFF_VDW_ENABLED] = false
 #   mm.options[:distance_dependent_dielectric] = true
#
 #   setup!(mm)
 #   nbc = mm.components[3]
 #   es = first(nbc.es_interactions)
 #   es.qi = -1
 #   es.qj = 2
#
 #   a2.r = Vector3{Float32}(11, 0, 0)
 #   atoms_df(syst).F .= Ref(Vector3(0,0,0))
 #   update!(nbc)
 #   es_ball = compute_energy(nbc)
 #   compute_forces(nbc)
#syst = load_sdfile(ball_data_path("../test/data/sdfile_test_1.sdf"))
#syst = deserialize("zsys")
#mm = MMFF94FF(syst)
#compute_energy(mm.components[3])
#compute_forces(mm.components[3])

#test_mm()

#sys = load_sdfile(ball_data_path("../test/data/MMFF94-stretch.sdf"))
#mm = MMFF94FF(sys)
#compute_forces(mm)
#e =compute_energy(mm)
#e2 = compute_energy(mm.components[1])

#test()






