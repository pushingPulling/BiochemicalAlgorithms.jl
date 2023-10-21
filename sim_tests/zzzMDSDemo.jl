using BiochemicalAlgorithms
BA = BiochemicalAlgorithms
using Plots


function verlet_at(syst,mm, len, step, dampen)   
    trajectories = []
    update!(mm.components[1])
    update!(mm.components[2])
    update!(mm.components[4])
    compute_forces(mm)
    nats = natoms(syst)
    forces = []
    force_now = Vector3(0,0,0)
    force_before = Vector3(0,0,0)
    force_after = Vector3(0,0,0)

    x_nm1 = []
    x_n = []
    x_np1 = Any[0 for _ in 1:nats]

    for (idx, at) in enumerate(atoms(syst))
        x_0 = at.r

        push!(x_nm1, x_0)
        push!(x_n, x_nm1[idx] + at.v * step + (1/2) * at.F * step^2)
        at.r = x_n[idx]

        push!(trajectories, [x_0, x_n[idx]])
    end

    at1 = atoms(syst)[1]
    force_now = at1.F


    for i in 1:len-2
        update!(mm.components[1])
        update!(mm.components[2])
        update!(mm.components[4])
        compute_forces(mm)
        for (idx, at) in enumerate(atoms(syst))
            x_np1[idx] = 2*x_n[idx] - x_nm1[idx] + at.F * step^2
            x_nm1[idx] = x_n[idx]
            x_n[idx] = x_np1[idx]
            at.r = x_n[idx]
            push!(trajectories[idx], x_n[idx])
        end
        force_after = at1.F
        abs(force_after.x) < abs(force_now.x) && abs(force_before.x) < abs(force_now.x) && (push!(forces, (i, at1.F.x)))
        force_before = force_now
        force_now = force_after

        i%(len/100) == 0 && println("$(100*i/len)% done")
    end

    foreach(println, forces)
    return trajectories
end

function verlet_const(step, len)
    A() = Vector3(0,0,-9.81)
    trajectory = []

    x_0 = Vector3(0,0,10)
    v_0 = Vector3(0,0,0)

    x_nm1 = x_0
    x_n = x_nm1 + v_0 * step + (1/2) * A() * step^2

    push!(trajectory, x_0, x_n)

    for _ in 1:len-2
        x_np1 = 2*x_n - x_nm1 + A() * step^2
        x_nm1 = x_n
        x_n = x_np1

        push!(trajectory, x_n)
    end

    return trajectory
end

#=
len = 31
x = verlet_const(0.05, len)
@gif for i in 1:len
    scatter([x[i][1]], [x[i][2]], [x[i][3]],  zlims = (0,20))
    #scatter!([-x[i][1]], [-x[i][2]], [20 - x[i][3]],  zlims = (0,20))
end
=#
syst = load_sdfile(ball_data_path("../test/data/MMFF94-stretch.sdf"))
mm = MMFF94FF(syst)
len = 1000
x = verlet_at(syst,mm, len, 10, 3.0)

@gif for i in 1:len
    plot([b[i].x for b in x[[1,2]]], [b[i].y for b in x[[1,2]]], [b[i].z for b in x[[1,2]]],
    marker=(:circle,5), xlim=(-1,2))
end every 10
@show x[1][end]
@show x[2][end];nothing