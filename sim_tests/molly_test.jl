using Molly
using Enzyme
using Zygote
using Flux
using Suppressor

function sim_random()
    n_atoms = 100
    boundary = CubicBoundary(2.0u"nm")
    temp = 298.0u"K"
    atom_mass = 10.0u"g/mol"

    atoms = [Atom(mass=atom_mass, σ=0.3u"nm", ϵ=0.2u"kJ * mol^-1") for i in 1:n_atoms]
    coords = place_atoms(n_atoms, boundary; min_dist=0.3u"nm")
    velocities = [random_velocity(atom_mass, temp) for i in 1:n_atoms]
    pairwise_inters = (LennardJones(),)
    simulator = VelocityVerlet(
        dt=0.002u"ps",
        coupling=AndersenThermostat(temp, 1.0u"ps"),
    )

    sys = System(
        atoms=atoms,
        coords=coords,
        boundary=boundary,
        velocities=velocities,
        pairwise_inters=pairwise_inters,
        loggers=(temp=TemperatureLogger(100),),
    )

    simulate!(sys, simulator, 10_000, n_threads = 1)
end

function sim_prot()

    sys = System(
        joinpath(dirname(pathof(Molly)), "..", "data", "5XER", "gmx_coords.gro"),
        joinpath(dirname(pathof(Molly)), "..", "data", "5XER", "gmx_top_ff.top");
        loggers=(
            temp=TemperatureLogger(10),
            writer=StructureWriter(10, "traj_5XER_1ps.pdb"),
        ),
    )

    temp = 298.0u"K"
    random_velocities!(sys, temp)
    simulator = VelocityVerlet(
        dt=0.0002u"ps",
        coupling=AndersenThermostat(temp, 1.0u"ps"),
    )

    simulate!(sys, simulator, 5_000,n_threads=1)
end

#sim_prot()


# function loss(sys_clean_energy, bond_k, bond_r0, neighbors, ff_dirty)
#     sys_dirty = System(
#         joinpath(data_dir, "6mrr_equil.pdb"),
#         ff_dirty;
#         velocities=sys_clean.velocities,
#         loggers=(
#             energy=TotalEnergyLogger(10),
#         ),
#         gpu=false,
#         units=false
#     )

#     simulate!(sys_dirty, simulator, 10; n_threads=1, run_loggers=false)
#     # return energy of system
#     return abs(total_energy(sys_dirty, neighbors) - sys_clean_energy) 
# end

#ToDo:
# use AD to get loss for ML
# have an unchanged forcefield for energy calc. Use dirty FF for sim and change it after each iter (ML)
# evaluate: How close is the new value to the original value?


# in the dirty version the "H/N3" Harmonic bond force is changed length 0.101->0.16; k 363171.19999999995->333171.19999999995
# the file is in C:\Users\Dan\.julia\packages\Molly\3eSlU\data\force_fields, line 4993

# get energy of unaltered simulation
# define simulation with slightly changed parameters for one bond.
# use AD to define loss function which compares energies of original and altered system
# use loss function for machine learning

# how di make Molly not use unitful?
# how do i apply the grad to the values?

# sys_clean energy = 61640.912760568724
# modified loss energy = 244440.0435986874          ~4 times higher
function AD_test_no_unit()
    data_dir = joinpath(dirname(pathof(Molly)), "..", "data") 
    ff_clean = MolecularForceField(
        joinpath(data_dir, "force_fields", "ff99SBildn.xml"),
        joinpath(data_dir, "force_fields", "tip3p_standard.xml"),
        joinpath(data_dir, "force_fields", "his.xml");
        units=false
    )

    sys_clean = System(
        joinpath(data_dir, "6mrr_equil.pdb"),
        ff_clean;
        units=false,
        gpu=false,
    )

    temp = 298.0
    simulator = Langevin(
        dt=0.001,    #u"ps",
        temperature=temp,
        friction=1.0,    #u"ps^-1",
        coupling=MonteCarloBarostat(1.0, temp, sys_clean.boundary),   #1.0 u"bar"
    )

    # this bond object holds the 2 parameters of interest
    bond_params  = sys_clean.specific_inter_lists[1].inters[1]
    bond_k_orig  = bond_params.k        # k  = 363171.19999999995
    bond_r0_orig = bond_params.r0       # r0 = 0.101

    random_velocities!(sys_clean, temp)

    sys_dirty      = deepcopy(sys_clean)
    sys_dirty_copy = deepcopy(sys_clean)

    #get energy from a run using unaltered parameters
    neighbors = find_neighbors(sys_clean, sys_clean.neighbor_finder; n_threads=8)
    simulate!(sys_clean, simulator, 10)
    sys_clean_energy = total_energy(sys_clean, neighbors)       # energy = 8632.53

    bond_k_alt  = Float64(100)
    bond_r0_alt = Float64(10)

    # run with unaltered parameters. Expect ≈0 loss and ≈0 derivative
    # returns ((0.0, 0.0, nothing, nothing, nothing),)
    grads_orig = autodiff(Reverse, loss, Active, Active(bond_k_orig), Active(bond_r0_orig), Const(sys_dirty), Const(sys_clean_energy), Const(simulator), Const(neighbors))
    sys_dirty = deepcopy(sys_dirty_copy)
    loss_orig = loss(bond_k_orig, bond_r0_orig, sys_dirty, sys_clean_energy, simulator, neighbors) # loss = 12.44
    sys_dirty = deepcopy(sys_dirty_copy)


    # run with   altered parameters. Expect high loss and high derivative
    # returns ((0.0, 0.0, nothing, nothing, nothing),)
    grads_alt = autodiff(Reverse, loss, Active, Active(bond_k_alt), Active(bond_r0_alt),   Const(sys_dirty), Const(sys_clean_energy), Const(simulator), Const(neighbors))
    sys_dirty = deepcopy(sys_dirty_copy)
    loss_alt  = loss(bond_k_alt, bond_r0_alt, sys_dirty, sys_clean_energy, simulator, neighbors)  # loss = 4907.96
    
    return grads_orig, loss_orig, grads_alt, loss_alt
end

function AD_test_enzyme()
    data_dir = joinpath(dirname(pathof(Molly)), "..", "data") 
    ff_clean = MolecularForceField(
        joinpath(data_dir, "force_fields", "ff99SBildn.xml"),
        joinpath(data_dir, "force_fields", "tip3p_standard.xml"),
        joinpath(data_dir, "force_fields", "his.xml");
        units=false
    )

    sys_clean = System(
        joinpath(data_dir, "6mrr_nowater.pdb"),
        ff_clean;
        units=false,
        gpu=false,
    )

    temp = 298.0
    simulator = Langevin(
        dt=0.001,    #u"ps",
        temperature=temp,
        friction=1.0,    #u"ps^-1",
        coupling=MonteCarloBarostat(1.0, temp, sys_clean.boundary),   #1.0 u"bar"
    )

    # this bond object holds the 2 parameters of interest
    bond_params  = sys_clean.specific_inter_lists[1].inters[1]
    bond_k_orig  = bond_params.k        # k  = 363171.19999999995
    bond_r0_orig = bond_params.r0       # r0 = 0.101

    random_velocities!(sys_clean, temp)

    sys_dirty      = deepcopy(sys_clean)
    sys_dirty_copy = deepcopy(sys_clean)

    #get energy from a run using unaltered parameters
    neighbors = find_neighbors(sys_clean, sys_clean.neighbor_finder; n_threads=8)
    simulate!(sys_clean, simulator, 10)
    sys_clean_energy = total_energy(sys_clean, neighbors)       # energy = 8632.53

    bond_k_alt  = Float64(100)
    bond_r0_alt = Float64(10)

    # run with unaltered parameters. Expect ≈0 loss and ≈0 derivative
    # returns ((0.0, 0.0, nothing, nothing, nothing),)
    grads_orig = autodiff(Reverse, loss, Active, Active(bond_k_orig), Active(bond_r0_orig), Active(sys_dirty), Const(sys_clean_energy), Const(simulator), Const(neighbors))
    sys_dirty = deepcopy(sys_dirty_copy)
    loss_orig = loss(bond_k_orig, bond_r0_orig, sys_dirty, sys_clean_energy, simulator, neighbors) # loss = 12.44
    sys_dirty = deepcopy(sys_dirty_copy)


    # run with   altered parameters. Expect high loss and high derivative
    # returns ((0.0, 0.0, nothing, nothing, nothing),)
    grads_alt = autodiff(Reverse, loss_enzyme, Active, Active(bond_k_alt), Active(bond_r0_alt),   Active(sys_dirty), Const(sys_clean_energy), Const(simulator), Const(neighbors))
    sys_dirty = deepcopy(sys_dirty_copy)
    loss_alt  = loss(bond_k_alt, bond_r0_alt, sys_dirty, sys_clean_energy, simulator, neighbors)  # loss = 4907.96
    
    return grads_orig, loss_orig, grads_alt, loss_alt
end

function loss_enzyme(k, r0, sys_dirty, sys_clean_energy, simulator, neighbors)
    sys_dirty.specific_inter_lists[1].inters[1] = HarmonicBond{Float64, Float64}(k, r0)
    simulate!(sys_dirty, simulator, 10; n_threads=1, run_loggers=false)
    return abs(total_energy(sys_dirty, neighbors) - sys_clean_energy)  
end


function AD_SD_Adam_zygote()
    data_dir = joinpath(dirname(pathof(Molly)), "..", "data") 
    ff_clean = MolecularForceField(
        joinpath(data_dir, "force_fields", "ff99SBildn.xml"),
        joinpath(data_dir, "force_fields", "tip3p_standard.xml"),
        joinpath(data_dir, "force_fields", "his.xml");
        units=false
    )

    sys_clean = System(
        joinpath(data_dir, "6mrr_equil.pdb"),
        ff_clean;
        units=false,
        gpu=false,
    )

    temp = 298.0
    simulator = Langevin(
        dt=0.001,    #u"ps",
        temperature=temp,
        friction=1.0,    #u"ps^-1",
        coupling=MonteCarloBarostat(1.0, temp, sys_clean.boundary),   #1.0 u"bar"
    )

    # this bond object holds the 2 parameters of interest
    bond_params  = sys_clean.specific_inter_lists[1].inters[1]
    bond_k_orig  = bond_params.k        # k  = 363171.19999999995
    bond_r0_orig = bond_params.r0       # r0 = 0.101

    random_velocities!(sys_clean, temp)

    sys_dirty      = deepcopy(sys_clean)

    neighbors = find_neighbors(sys_clean, sys_clean.neighbor_finder; n_threads=8)
    simulate!(sys_clean, simulator, 10)
    sys_clean_energy = total_energy(sys_clean, neighbors)       # energy = 8632.53

    bond_k_alt  = Float64(500) # bond_k_orig
    bond_r0_alt = Float64(3)   # bond_r0_orig

    x = [bond_k_alt, bond_r0_alt]
    opt =  Flux.Optimisers.Adam(0.9, (0.9, 0.999), 1.0e-8)
    leaf = Flux.Optimisers.Leaf(opt, (zeros(2),zeros(2),opt.beta), false)

    n_epochs = 100

    

    for epoch_n in 1:n_epochs
        grad_k, grad_r0, _, _, _, _ = Zygote.gradient(loss, bond_k_alt, bond_r0_alt, sys_dirty,sys_clean_energy, simulator, neighbors)
        loss_cur = loss(bond_k_alt, bond_r0_alt, sys_dirty,sys_clean_energy, simulator, neighbors)

        st, dx = Flux.Optimisers.apply!(leaf.rule, leaf.state, Float64[], [grad_k, grad_r0])
        leaf.state = st
        Flux.Optimisers.subtract!(x, dx)
        
        println("Epoch ", epoch_n, " loss ",loss_cur ," k:", x[1], " r0:",x[2], " grad k:", grad_k, " grad r0: ", grad_r0)
    end
    # TodO: 
    # use optimizer object to optimize the gradients
    # Adam holds its poarams, its state.
    # as dx use the grad dl/dw: Meaning, for each param i wanna learn, use one leaf of the opt
    # train teh weight as such: w_t+1 = w_t - alpha*m_t
    # Todo: Figure out the leaf system
    # figure out how to use the params dict
    # unrelated: Actually try out a deep neural net with MDS
end


function AD_SD_Adam_zygote_all_bonds()
    data_dir = joinpath(dirname(pathof(Molly)), "..", "data") 
    ff_clean = MolecularForceField(
        joinpath(data_dir, "force_fields", "ff99SBildn.xml"),
        joinpath(data_dir, "force_fields", "tip3p_standard.xml"),
        joinpath(data_dir, "force_fields", "his.xml");
        units=false
    )

    sys_clean = System(
        joinpath(data_dir, "6mrr_equil.pdb"),
        ff_clean;
        units=false,
        gpu=false,
    )

    temp = 298.0
    simulator = Langevin(
        dt=0.001,    #u"ps",
        temperature=temp,
        friction=1.0,    #u"ps^-1",
        coupling=MonteCarloBarostat(1.0, temp, sys_clean.boundary),   #1.0 u"bar"
    )

    # this bond object holds the 2 parameters of interest
    bond_params = sys_clean.specific_inter_lists[1].inters
    #bond_params = vcat(collect([p.k,p.r0] for p in bond_params)...)

    # these random numbers may cause Floating point inaccuracy (NaN) issues; maybe changing the bounds will make this more stable
    bounds_k = 250_000.0:10:580_000.0
    bounds_r = 0.1:0.001:1.65
    params = vcat(map((x,y)-> [x,y], rand(bounds_k, length(bond_params)), rand(bounds_r, length(bond_params)))...)

    random_velocities!(sys_clean, temp)

    sys_dirty      = deepcopy(sys_clean)

    neighbors = find_neighbors(sys_clean, sys_clean.neighbor_finder; n_threads=8)
    simulate!(sys_clean, simulator, 10)
    sys_clean_energy = total_energy(sys_clean, neighbors)       # energy = 8632.53


    x = params
    opt =  Flux.Optimisers.Adam(0.02, (0.9, 0.999), 1.0e-8)
    leaf = Flux.Optimisers.Leaf(opt, (zeros(length(params)),zeros(length(params)),opt.beta), false)

    n_epochs = 10

    for epoch_n in 1:n_epochs
        grads, _, _, _, _ = Zygote.jacobian(loss_all_bonds, params, sys_dirty,sys_clean_energy, simulator, neighbors)
        loss_cur = loss_all_bonds( params, sys_dirty,sys_clean_energy, simulator, neighbors)

        st, dx = Flux.Optimisers.apply!(leaf.rule, leaf.state, Float64[], grads' )
        leaf.state = st
        Flux.Optimisers.subtract!(x, dx)
        
        println("Epoch ", epoch_n, " loss ",loss_cur ," k:", x[1], " r0:",x[2])
    end
    # TodO: 
    # use optimizer object to optimize the gradients
    # Adam holds its poarams, its state.
    # as dx use the grad dl/dw: Meaning, for each param i wanna learn, use one leaf of the opt
    # train teh weight as such: w_t+1 = w_t - alpha*m_t
    # Todo: Figure out the leaf system
    # figure out how to use the params dict
    # unrelated: Actually try out a deep neural net with MDS
end

function AD_test_zygote()
    data_dir = joinpath(dirname(pathof(Molly)), "..", "data") 
    ff_clean = MolecularForceField(
        joinpath(data_dir, "force_fields", "ff99SBildn.xml"),
        joinpath(data_dir, "force_fields", "tip3p_standard.xml"),
        joinpath(data_dir, "force_fields", "his.xml");
        units=false
    )

    sys_clean = System(
        joinpath(data_dir, "6mrr_nowater.pdb"),
        ff_clean;
        units=false,
        gpu=false,
    )

    temp = 298.0
    simulator = Langevin(
        dt=0.001,    #u"ps",
        temperature=temp,
        friction=1.0,    #u"ps^-1",
        coupling=MonteCarloBarostat(1.0, temp, sys_clean.boundary),   #1.0 u"bar"
    )

    # this bond object holds the 2 parameters of interest
    bond_params  = sys_clean.specific_inter_lists[1].inters[1]
    bond_k_orig  = bond_params.k        # k  = 363171.19999999995
    bond_r0_orig = bond_params.r0       # r0 = 0.101

    random_velocities!(sys_clean, temp)

    sys_dirty      = deepcopy(sys_clean)
    sys_dirty_copy = deepcopy(sys_clean)

    #get energy from a run using unaltered parameters
    neighbors = find_neighbors(sys_clean, sys_clean.neighbor_finder; n_threads=8)
    simulate!(sys_clean, simulator, 10)
    sys_clean_energy = total_energy(sys_clean, neighbors)       # energy = 8632.53

    #alter the parameters
    bond_k_alt  = bond_k_orig * 0.6
    bond_r0_alt = bond_r0_orig * 0.6

    σlearn_k = 1000.00000
    σlearn_r0 = 0.000010
    n_epochs = 200
    println("bond_k $(bond_k_alt) bond_r0 $(bond_r0_alt)")

    for epoch_n in 1:n_epochs
        grad_k, grad_r0, _, _, _, _ = Zygote.gradient(loss_copy_and_buf, bond_k_alt, bond_r0_alt, sys_dirty,sys_clean_energy, simulator, neighbors)
        sys_dirty = deepcopy(sys_dirty_copy)
        loss=loss_copy_and_buf(bond_k_alt, bond_r0_alt, sys_dirty,sys_clean_energy, simulator, neighbors)
        sys_dirty = deepcopy(sys_dirty_copy)

        bond_k_alt -= grad_k * 1e-2 * σlearn_k
        bond_r0_alt -= grad_r0 * 7e-2 * σlearn_r0
        println("Epoch ", epoch_n, " loss ",loss ," k:", bond_k_alt, " r0:",bond_r0_alt, " grad k:", grad_k, " grad r0: ", grad_r0)
    end
    
    return bond_k_alt, bond_r0_alt
end


function loss(k, r0, sys_dirty_orig, sys_clean_energy, simulator, neighbors)
    sys_dirty2 = Zygote.ignore() do 
        deepcopy(sys_dirty_orig)
    end
    inner = InteractionList2Atoms{Vector{Int32}, Vector{HarmonicBond{Float64, Float64}}}( sys_dirty2.specific_inter_lists[1].is, sys_dirty2.specific_inter_lists[1].js, vcat(HarmonicBond{Float64, Float64}(k, r0), sys_dirty2.specific_inter_lists[1].inters[2:end]), sys_dirty2.specific_inter_lists[1].types)
    lis = (inner, sys_dirty2.specific_inter_lists[2:end]...)
    sys_dirty2.specific_inter_lists = lis
    simulate!(sys_dirty2, simulator, 10; n_threads=1, run_loggers=false)
    return abs(total_energy(sys_dirty2, neighbors, n_threads=1) - sys_clean_energy)  
end


function loss_jac(params, sys_dirty_orig, sys_clean_energy, simulator, neighbors)
    sys_dirty2 = Zygote.ignore() do 
        deepcopy(sys_dirty_orig)
    end
    k, r0 = params
    inner = InteractionList2Atoms{Vector{Int32}, Vector{HarmonicBond{Float64, Float64}}}( sys_dirty2.specific_inter_lists[1].is, sys_dirty2.specific_inter_lists[1].js, vcat(HarmonicBond{Float64, Float64}(k, r0), sys_dirty2.specific_inter_lists[1].inters[2:end]), sys_dirty2.specific_inter_lists[1].types)
    lis = (inner, sys_dirty2.specific_inter_lists[2:end]...)
    sys_dirty2.specific_inter_lists = lis
    simulate!(sys_dirty2, simulator, 10; n_threads=1, run_loggers=false)
    return abs(total_energy(sys_dirty2, neighbors, n_threads=1) - sys_clean_energy)  
end

function loss_all_bonds(params, sys_dirty_orig, sys_clean_energy, simulator, neighbors)
    sys_dirty2 = Zygote.ignore() do 
        deepcopy(sys_dirty_orig)
    end
    all_bonds = vcat(collect(HarmonicBond{Float64,Float64}(params[i], params[i+1]) for i in 1:2:length(params)-1))
    inner = InteractionList2Atoms{Vector{Int32}, Vector{HarmonicBond{Float64, Float64}}}( sys_dirty2.specific_inter_lists[1].is, sys_dirty2.specific_inter_lists[1].js, all_bonds, sys_dirty2.specific_inter_lists[1].types)
    lis = (inner, sys_dirty2.specific_inter_lists[2:end]...)
    sys_dirty2.specific_inter_lists = lis
    simulate!(sys_dirty2, simulator, 10; n_threads=1, run_loggers=false)
    return abs(total_energy(sys_dirty2, neighbors, n_threads=1) - sys_clean_energy)  
end


function loss_copy_and_buf(k, r0, sys_dirty, sys_clean_energy, simulator, neighbors)
    sys_dirty2 = deepcopy(sys_dirty);    inner = InteractionList2Atoms{Vector{Int32}, Vector{HarmonicBond{Float64, Float64}}}( sys_dirty2.specific_inter_lists[1].is, sys_dirty2.specific_inter_lists[1].js, vcat(HarmonicBond{Float64, Float64}(k, r0), sys_dirty2.specific_inter_lists[1].inters[2:end]), sys_dirty2.specific_inter_lists[1].types);    lis = (inner, sys_dirty2.specific_inter_lists[2:end]...);    sys_dirty2.specific_inter_lists = lis;
    simulate!(sys_dirty2, simulator, 10; n_threads=1, run_loggers=false)
    return abs(total_energy(sys_dirty2, neighbors, n_threads=1) - sys_clean_energy)  
end




function AD_test(starting_learning_rate)

    # ----------------------------------
    #define sim and loss func which will be enyzme'd
    data_dir = joinpath(dirname(pathof(Molly)), "..", "data") 
    ff_clean = MolecularForceField(
        joinpath(data_dir, "force_fields", "ff99SBildn.xml"),
        joinpath(data_dir, "force_fields", "tip3p_standard.xml"),
        joinpath(data_dir, "force_fields", "his.xml"),
    )

    ff_dirty = MolecularForceField(
        joinpath(data_dir, "force_fields", "ff99SBildn.xml"),
        joinpath(data_dir, "force_fields", "tip3p_standard.xml"),
        joinpath(data_dir, "force_fields", "his.xml"),
    )
    #alter ff - change parameters by ~58%
    bond_type_h_n3      = ff_dirty.bond_types[("H","N3")]
    k_type = typeof(bond_type_h_n3.k)        #original val: 363171.19999999995
    r0_type = typeof(bond_type_h_n3.r0)         #original val 0.101
    bond_type_h_n3 = typeof(bond_type_h_n3)(k_type(229251.81999999998), r0_type(0.160))


    
    sys_clean = System(
        joinpath(data_dir, "6mrr_nowater.pdb"),
        ff_clean;
        loggers=(
            energy=TotalEnergyLogger(10),
        ),
        gpu=false,
    )

    temp = 298.0u"K"
    random_velocities!(sys_clean, temp)

    sys_dirty = System(
        joinpath(data_dir, "6mrr_nowater.pdb"),
        ff_clean;
        velocities=sys_clean.velocities,
        loggers=(
            energy=TotalEnergyLogger(10),
        ),
        gpu=false,
    )

    simulator = Langevin(
        dt=0.001u"ps",    #u"ps",
        temperature=temp,
        friction=1.0u"ps^-1",    #u"ps^-1",
        coupling=MonteCarloBarostat(1.0u"bar", temp, sys_clean.boundary),   #1.0 u"bar"
    )
    neighbors = find_neighbors(sys_clean, sys_clean.neighbor_finder; n_threads=8)
    sys_clean_energy = total_energy(sys_clean, neighbors)
    simulate!(sys_clean, simulator, 500)
    sys_clean_energy = total_energy(sys_clean, neighbors)

    #loss func

    function loss(sys_dirty, sys_clean_energy, bond_k, bond_r0)
        bond_type_h_n3      = ff_dirty.bond_types[("H","N3")]
        k_type = typeof(bond_type_h_n3.k)        #original val: 363171.19999999995
        r0_type = typeof(bond_type_h_n3.r0)         #original val 0.101
        bond_type_h_n3 = typeof(bond_type_h_n3)(k_type(bond_k), r0_type(bond_r0))

        simulate!(sys, simulator, 500)
        # return energy of system
        return abs(total_energy(sys_dirty) - sys_clean_energy) 
    end

    bond_k  = ustrip(bond_type_h_n3.k )
    bond_r0 = ustrip(bond_type_h_n3.r0)
    σlearn = starting_learning_rate
    n_epochs = 10

    for epoch_n in 1:n_epochs
        #grad = gradient(loss, σlearn, coords, velocities)[1]
        #call loss here
        #  |  |  ||
        # ____|____
        #     |
        #  || |  | _
        #reverse accu https://enzyme.mit.edu/index.fcgi/julia/stable/
        # seems fine up until here. Need to find a way to make Molly go completely without Unitful
        grads = autodiff(Reverse, loss, Active, Const(sys_dirty), Const(sys_clean_energy), Active(bond_k), Active(bond_r0))
        println("Epoch $(epoch_n) | Grad $(round(grad,digits=5))\n", grad)
        σlearn -= abs(grad * 1e-2)
    end
    

    
end

function f(x,y,z)
    return x*x + y*y + z *z
end

function f2(x,y,z)
    return x*x + y*y + z *z, x^3+y^3+z^3
end

function f2_rev(x,y,z, res)
    push!(res,x*x + y*y + z *z)
    push!(res, x^3+y^3+z^3)
    return nothing
end

function test_enzyme()
    # f R^n -> R^m
    # possibly wrong explanations ahead
    # see file:///C:/Users/Dan/Desktop/uni/papers%20for%20Master%20thesis/Differentiable%20molecular%20simulation%20can%20learn.pdf
    #    for short description  on page 11

    # ?
    # use reverse for n >> m. useful to compute df/dx for each x that is wished
    # only one return, or return nothign and store output in a list passed as argmunet
    grad = autodiff(Reverse, f, Active, Active(2), Active(3), Active(4))
    @show grad
    #reverse for multiple returns: put all results into a list and return nothing
    res = []
    grad = autodiff(Reverse, f2_rev, Const, Active(2), Active(3), Active(4), Const(res))
    @show grad
    @show res

    # ?
    #use forward for n << m. useful to compute sum of all df/dx for each x that is wished (ans also if f has multiple returns)
    # need to duplicated because im carryin g the values to the end (and then reverse pass?)
    # reutrn is ((normal run returns) , (derivative returns))
    grad = autodiff(Forward, f, Duplicated, Duplicated(2.0,1.0),Duplicated(3.0,1.0), Const(4.0))
    @show grad

end

function time_vanilla()


end