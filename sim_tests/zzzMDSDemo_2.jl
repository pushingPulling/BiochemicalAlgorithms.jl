using BiochemicalAlgorithms
BA = BiochemicalAlgorithms
using Plots
using PlotlyJS
using Statistics
using Serialization
using LinearAlgebra

sys = load_sdfile("$(pwd())/test/data/rings_test.sdf")
mm = MMFF94FF(sys; tree=true)


function verlet_at(syst,mm, len, step, comps)

    trajectories = []
    BA.update!(mm)
    compute_forces(mm)
    nats = natoms(syst)
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

    for i in 1:len-2
        BA.update!(mm)
        compute_forces(mm)
        for (idx, at) in enumerate(atoms(syst))
            x_np1[idx] = 2*x_n[idx] - x_nm1[idx] + at.F * step^2
            x_nm1[idx] = x_n[idx]
            x_n[idx] = x_np1[idx]
            at.r = x_n[idx]
            push!(trajectories[idx], x_n[idx])
        end
        #i%(len/100) == 0 && println("$(100*i/len)% done")
    end
    return trajectories
end

function verlet_at_stat(syst,mm, len, step, comps, free)

    at_orig = deepcopy(atoms_df(syst).r)
    @show atoms_df(syst).number
    trajectories = []
    BA.update!(mm)
    compute_forces(mm)
    nats = natoms(syst)
    x_nm1 = []
    x_n = []
    x_np1 = Any[0 for _ in 1:nats]
    for (idx, at) in enumerate(atoms(syst))
        x_0 = at.r
        push!(x_nm1, x_0)
        push!(x_n, x_nm1[idx] + at.v * step + (1/2) * at.F * step^2)
        at.r = x_n[idx]
        if at.number == free
            @show at.number
            at.r = x_n[idx] =  at_orig[at.number]
        end
        push!(trajectories, [x_0, x_n[idx]])
    end

    for i in 1:len-2
        BA.update!(mm)
        compute_forces(mm)
        for (idx, at) in enumerate(atoms(syst))
            x_np1[idx] = 2*x_n[idx] - x_nm1[idx] + at.F * step^2
            x_nm1[idx] = x_n[idx]
            x_n[idx] = x_np1[idx]
            at.r = x_n[idx]
            if at.number == free
                at.r = x_n[idx] = at_orig[at.number]
            end
            push!(trajectories[idx], x_n[idx])
        end
        #i%(len/100) == 0 && println("$(100*i/len)% done")
    end
    return trajectories
end

function verlet_at_comps(syst,mm, len, step, comps)
    a = atoms(syst)[2]
    ao = deepcopy(atoms(syst)[2].r) 
    trajectories = []
    map(BA.update!, mm.components[comps])
    map(compute_forces, mm.components[comps])
    #BA.update!(mm)
    #compute_forces(mm)
    nats = natoms(syst)
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

    for i in 1:len-2
        map(BA.update!, mm.components[comps])
        map(compute_forces, mm.components[comps])
        #BA.update!(mm)
        #compute_forces(mm)
        for (idx, at) in enumerate(atoms(syst))

            x_np1[idx] = 2*x_n[idx] - x_nm1[idx] + at.F * step^2
            x_nm1[idx] = x_n[idx]
            x_n[idx] = x_np1[idx]
            # if a.number == at.number
            #     x_n[idx] = ao
            # end
            at.r = x_n[idx]
            #a.r = ao
            push!(trajectories[idx], x_n[idx])
        end
        i%(len/100) == 0 && println("$(100*i/len)% done")
    end
    return trajectories
end

function verlet_at_bend(syst,mm, len, step, comps)
    at1,at2,at3 = atoms(syst)
    a(idx) = x_n[idx]
    d21 = norm(at2.r-at1.r)
    d23 = norm(at2.r-at3.r)
    a = atoms(syst)[2]
    ao = deepcopy(atoms(syst)[2].r) 
    trajectories = []
    # map(update!, mm.components[comps])
    # map(compute_forces, mm.components[comps])
    BA.update!(mm)
    compute_forces(mm)
    nats = natoms(syst)
    x_nm1 = []
    x_n = []
    x_np1 = Any[0 for _ in 1:nats]
    for (idx, at) in enumerate(atoms(syst))
        x_0 = at.r
        push!(x_nm1, x_0)
        push!(x_n, x_nm1[idx] + at.v * step + (1/2) * at.F * step^2)
        if 1 == at.number
            #x_n[1] = at2.r + (at2.r - x_n[1]) * d21/norm(at2.r - x_n[1])
           # at1.r = x_n[1]
        end 
        if 3 == at.number
            #x_n[3] = at2.r + (at2.r - x_n[3]) * d23/norm(at2.r - x_n[3])
            #at3.r = x_n[3]
        end 
        if a.number == at.number
            #x_n[2] = Vector3(0,0,0)
            #at.r = Vector3(0,0,0)
        end
        
        #a.r = ao
        push!(trajectories, [x_0, x_n[idx]])
    end

    for i in 1:len-2
        # map(update!, mm.components[comps])
        # map(compute_forces, mm.components[comps])
        BA.update!(mm)
        compute_forces(mm)
        for (idx, at) in enumerate(atoms(syst))
            x_np1[idx] = 2*x_n[idx] - x_nm1[idx] + at.F * step^2
            x_nm1[idx] = x_n[idx]
            x_n[idx] = x_np1[idx]
            at.r = x_n[idx]
            if a.number == at.number
                #x_n[2] = at.r = Vector3(0,0,0)
            end
            if 1 == at.number
                #x_n[1] = x_n[1] * d21/norm(x_n[1])
                #at1.r = x_n[1]
            end 
            if 3 == at.number
                #x_n[3] = x_n[3] * d23/norm(x_n[3])
                #at3.r = x_n[3]
            end  
            push!(trajectories[idx], x_n[idx])
        end
        i%(len/100) == 0 && println("$(100*i/len)% done")
    end
    return trajectories
end

function verlet_at_strebend(syst,mm, len, step, comps)
    #a is the stem, b is being projected
    projection(a,b) = dot(a,b)*a/norm(a)^2
    at1,at2,at3 = atoms(syst)
    a(idx) = x_n[idx]
    v21 = at2.r-at1.r
    v23 = at2.r-at3.r
    a = atoms(syst)[2]
    ao = deepcopy(atoms(syst)[2].r) 
    trajectories = []
    # map(update!, mm.components[comps])
    # map(compute_forces, mm.components[comps])
    BA.update!(mm)
    compute_forces(mm)
    nats = natoms(syst)
    x_nm1 = []
    x_n = []
    x_np1 = Any[0 for _ in 1:nats]
    for (idx, at) in enumerate(atoms(syst))
        x_0 = at.r
        push!(x_nm1, x_0)
        push!(x_n, x_nm1[idx] + at.v * step + (1/2) * at.F * step^2)
        if 1 == at.number
            x_n[1] = projection(v21, x_n[1])
            at1.r = x_n[1]
        end 
        if 3 == at.number
            x_n[3] = projection(v23, x_n[3])
            at3.r = x_n[3]
        end 
        if a.number == at.number
            x_n[2] = Vector3(0,0,0)
            at.r = Vector3(0,0,0)
        end
        
        #a.r = ao
        push!(trajectories, [x_0, x_n[idx]])
    end

    for i in 1:len-2
        # map(update!, mm.components[comps])
        # map(compute_forces, mm.components[comps])
        BA.update!(mm)
        compute_forces(mm)
        for (idx, at) in enumerate(atoms(syst))
            x_np1[idx] = 2*x_n[idx] - x_nm1[idx] + at.F * step^2
            x_nm1[idx] = x_n[idx]
            x_n[idx] = x_np1[idx]
            at.r = x_n[idx]
            if a.number == at.number
                x_n[2] = at.r = Vector3(0,0,0)
            end
            if 1 == at.number
                x_n[1] = projection(v21, x_n[1])
                at1.r = x_n[1]
            end 
            if 3 == at.number
                x_n[3] = projection(v23, x_n[3])
                at3.r = x_n[3]
            end 
            push!(trajectories[idx], x_n[idx])
        end
        i%(len/100) == 0 && println("$(100*i/len)% done")
    end
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


function stretch(len, step, every_)
    syst = load_sdfile(ball_data_path("../test/data/MMFF94-stretch.sdf"))
    mm = MMFF94FF(syst)
    at1, at2 = atoms(syst)
    at2.r = Vector3(1.3720001,0,0)
    BA.update!(mm)
    x = verlet_at(syst, mm, len, step, [1,2,4])
    @show mean(x[1])
    @show mean(x[2])
#
    @gif for i in 1:len
        #scatter([x[1][i].x], [x[1][i].y], [x[1][i].z], xlim=(-2,2))
        #scatter!([x[2][i].x], [x[2][i].y], [x[2][i].z])
        Plots.plot([x[1][i].x, x[2][i].x], [x[1][i].y, x[2][i].y], [x[1][i].z, x[2][i].z], xlim=(-2,2),marker=(:circle,5), label=nothing, dpi=600)
    end every every_
end

#stretch(1000,50,10)

rotz(a) = [cos(deg2rad(a)) -sin(deg2rad(a)) 0; sin(deg2rad(a)) cos(deg2rad(a)) 0; 0 0 1]
function full_bend(len, step, every_, le, ro)
    syst = load_sdfile(ball_data_path("../test/data/MMFF94-bend.sdf"))
    mm = MMFF94FF(syst)
    mm.options[:stretches_enabled] = false
    mm.options[:bends_enabled] = true
    mm.options[:stretch_bend_enabled] = false
    at1, at2, at3 = atoms(syst)

    at3.r = rotz(ro) * at3.r
    @show rad2deg(acos(dot(at1.r,rotz(ro)*at3.r)/norm(at1.r)*(norm(at3.r))))
    
    v12 = at1.r-at2.r
    v23 = at3.r-at2.r
    at1.r = v12 * le
    at3.r = v23 * le
    BA.update!(mm)
    x = verlet_at_bend(syst, mm, len, step, [1,2,4])     #or verlet_at_bend
    @show mean(x[1])
    @show mean(x[2])
#
    @gif for i in 2:len
        Plots.plot([x[1][i].x, x[2][i].x], [x[1][i].y, x[2][i].y], [x[1][i].z, x[2][i].z],xlim=(-2,2), ylim=(-2,2), marker=(:circle,5), labels=nothing)
        Plots.plot!([x[3][i].x, x[2][i].x], [x[3][i].y, x[2][i].y], [x[3][i].z, x[2][i].z],xlim=(-2,2), ylim=(-2,2), marker=(:circle,5))
    end every every_
    #@show norm(at3.r - at1.r)
end


function bend(len, step, every_, ro)
    syst = load_sdfile(ball_data_path("../test/data/MMFF94-bend.sdf"))
    mm = MMFF94FF(syst)
    mm.options[:stretches_enabled] = true
    mm.options[:stretch_bends_enabled] = false
    at1, at2, at3 = atoms(syst)
    at1.r = rotz(-ro) * at1.r
    at3.r = rotz(ro) * at3.r
    BA.update!(mm)
    x = verlet_at_bend(syst, mm, len, step, [1,2,4])
    angle(a,b,g) = rad2deg(acos(dot(a-b,g-b)/(norm(a-b)*(norm(g-b)))))
    @show mean(x[1])
    @show mean(x[2])
    angles = map(angle, x...)
    @show mean(angles)
#
    @gif for i in 1:len
        # scatter([x[1][i].x], [x[1][i].y], [x[1][i].z], xlim=(-2,2), ylim=(-2,2))
        # scatter!([x[2][i].x], [x[2][i].y], [x[2][i].z])
        # scatter!([x[3][i].x], [x[3][i].y], [x[3][i].z])
        Plots.plot([x[1][i].x, x[2][i].x], [x[1][i].y, x[2][i].y], [x[1][i].z, x[2][i].z],xlim=(-2,2), ylim=(-2,2), marker=(:circle,5), legend=nothing, dpi=600)
        Plots.plot!([x[3][i].x, x[2][i].x], [x[3][i].y, x[2][i].y], [x[3][i].z, x[2][i].z],xlim=(-2,2), ylim=(-2,2), marker=(:circle,5))
    end every every_
end


#min angle -1.8808849831992736
#norm(rotz(-1.8808849831992736)*at3.r - at1.r) = 1.3926803635672056
# k = d/dot(a, a* norm(b)/norm(a) * cos(deg2rad({degree between a and b})))
# do norm(rotz(-30)*at3.r - at1.r) * f   and change f until this is = 1.3926
# values obtained like this also kinda work for verlet_at

#finally.
# full_bend(1000,50,10, 1.0,-14.90706) refeerence angle. change 
# full_bend(1000,50,10, 1.0,-14.90706) refeerence len. change angle to anything an nothing happens 
function strebend(len, step, every_, le, ro)
    syst = load_sdfile(ball_data_path("../test/data/MMFF94-bend.sdf"))
    mm = MMFF94FF(syst)
    mm.options[:stretches_enabled] = false       #this makes it look awesome
    mm.options[:bends_enabled] = false
    at1, at2, at3 = atoms(syst)
    at3.r = rotz(ro) * at3.r
    v12 = at1.r-at2.r
    v23 = at3.r-at2.r
    at1.r = v12 * le
    at3.r = v23 * le
    BA.update!(mm)
    x = verlet_at_strebend(syst, mm, len, step, [1,2,4])
    @show mean(x[1])
    @show mean(x[3])
#
    @gif for i in 1:len
        Plots.plot([x[1][i].x, x[2][i].x], [x[1][i].y, x[2][i].y], [x[1][i].z, x[2][i].z],xlim=(-4,4), ylim=(-4,4), marker=(:circle,5), labels=nothing, legend=nothing,dpi=600)
        Plots.plot!([x[3][i].x, x[2][i].x], [x[3][i].y, x[2][i].y], [x[3][i].z, x[2][i].z],  marker=(:circle,5))
    end every every_
    #return x
end

function plane(len, step, every_)
    syst = load_sdfile(ball_data_path("../test/data/MMFF94-plane.sdf"))
    at1,at2,at3,at4 = atoms(syst)
#
#
    mm = MMFF94FF(syst)
    @show atoms_df(syst).idx
    @show mm.components[4].out_of_plane_bends[1].i.idx
    @show mm.components[4].out_of_plane_bends[1].j.idx
    @show mm.components[4].out_of_plane_bends[1].k.idx
    @show mm.components[4].out_of_plane_bends[1].l.idx
    x = verlet_at(syst,mm, len, step, [4])
    ats = [at1.r, at3.r, at4.r]
#
    @gif for i in 1:len
        plot([b[i].x for b in x[[1,2]]], [b[i].y for b in x[[1,2]]], [b[i].z for b in x[[1,2]]],
            marker=(:circle,5), xlim=(-2,3), ylim=(-2,2), zlim=(-2,2))
        plot!([b[i][1] for b in x[2:3]], [b[i][2] for b in x[2:3]], [b[i][3] for b in x[2:3]],marker=(:circle,5))
        plot!([b[i][1] for b in x[[2,4]]], [b[i][2] for b in x[[2,4]]], [b[i][3] for b in x[[2,4]]],marker=(:circle,5))
        ats = [x[1][i], x[3][i], x[4][i]]
        m = mean(ats)
        n = cross(ats[2]- ats[1], ats[3] - ats[1]) * 0.3
    
        scatter!([m.x], [m.y], [m.z])
        quiver!([m.x], [m.y], [m.z], quiver=([n.x], [n.y], [n.z]))
        #plot!([m.x, (m-n).x], [m.y,(m-n).y], [m.z,(m-n).z], marker=(:circle,10))
        #
    end every every_
    return x[1][600], x[3][600], x[4][600], x[2][600]
end

#plane(1000,50,10)

function tors(len, step, every_)
    syst = load_sdfile(ball_data_path("../test/data/MMFF94-torsion.sdf"))
    at1,at2,at3,at4 = atoms(syst)
#
#
    mm = MMFF94FF(syst)
    x = verlet_at(syst,mm, len, step, [2])
#
    @gif for i in 1:len
        Plots.plot([b[i].x for b in x[[1,2]]], [b[i].y for b in x[[1,2]]], [b[i].z for b in x[[1,2]]],
            marker=(:circle,5), xlim=(-2,2), ylim=(-2,2), zlim=(-2,2))
        Plots.plot!([b[i][1] for b in x[2:3]], [b[i][2] for b in x[2:3]], [b[i][3] for b in x[2:3]],marker=(:circle,5))
        Plots.plot!([b[i][1] for b in x[[3,4]]], [b[i][2] for b in x[[3,4]]], [b[i][3] for b in x[[3,4]]],marker=(:circle,5))
        #
    end every every_
end
#tors(1000,18,10)

function vdw(len, step, every_)
    syst = load_sdfile(ball_data_path("../test/data/MMFF94-vdw.sdf"))
    at1,at2 = atoms(syst)#
    at2.r = Vector3(7,0,0)
#
#
    mm = MMFF94FF(syst)
    #mm.options[:MMFF_ES_ENABLED] = false
    x = verlet_at(syst,mm, len, step, [4])
#
    @gif for i in 1:len
        #scatter([b[i].x for b in x[1]], [b[i].y for b in x[1]], [b[i].z for b in x[1]],
         #   xlim=(-2,2), ylim=(-2,2), zlim=(-2,2))
        Plots.scatter([x[1][i].x], [x[1][i].y], [x[1][i].z], xlim=(-5,8),dpi=600, msize=5, color="red")
        Plots.scatter!([x[2][i].x], [x[2][i].y], [x[2][i].z], msize=5, color="red")
        #
        #
    end every every_
end
#vdw(1000,100,10)

function es(len, step, every_)
    syst = load_sdfile(ball_data_path("../test/data/MMFF94-vdw.sdf"))
    at1,at2 = atoms(syst)#
    at2.r = Vector3(5,0,0)
#
#
    mm = MMFF94FF(syst)
    mm.options[:MMFF_VDW_ENABLED] = false
    x = verlet_at(syst,mm, len, step, [4])
#
    @gif for i in 1:len
        #scatter([b[i].x for b in x[1]], [b[i].y for b in x[1]], [b[i].z for b in x[1]],
         #   xlim=(-2,2), ylim=(-2,2), zlim=(-2,2))
        Plots.scatter([x[1][i].x], [x[1][i].y], [x[1][i].z], xlim=(-6,7), color="red",msize=5, dpi=600)
        Plots.scatter!([x[2][i].x], [x[2][i].y], [x[2][i].z], color="red",msize=5)
        #
    end every every_
end

#es(5000,50,10)

function molec(reps, every_;len=100, name="Structure2D_COMPOUND_CID_196182.sdf")
    load_op = endswith(name,".sdf") ? load_sdfile : load_pdb
    syst = load_op(name, Float64)
    mm = MMFF94FF(syst)
    @show "done mm"
    j = [begin println("$(i)/$(reps) energy=$(compute_energy(mm))"); verlet_at(syst,mm,len,10,[]) end for i in 1:reps]
    x = [vcat((l[i] for l in j)...) for i in 1:natoms(syst)]
    
    xl = extrema(Iterators.flatten((coords.x for coords in  at) for at in  x))
    xl = (floor(xl[1]), ceil(xl[2]))
    yl = extrema(Iterators.flatten((coords.y for coords in  at) for at in  x))
    xl = (floor(yl[1]), ceil(yl[2]))
    zl = extrema(Iterators.flatten((coords.z for coords in  at) for at in  x))
    zl = (floor(zl[1]), ceil(zl[2]))

    @gif for i in 1:length(x[1])
        #scatter([b[i].x for b in x[1]], [b[i].y for b in x[1]], [b[i].z for b in x[1]],
         #   xlim=(-2,2), ylim=(-2,2), zlim=(-2,2))
        Plots.plot(zlim=zl, xlim=xl, ylim=yl)
        #plot()
        for bond in bonds(syst) 
            Plots.plot!([b[i].x for b in x[[bond.a1-1,bond.a2-1]]], [b[i].y for b in x[[bond.a1-1,bond.a2-1]]], [b[i].z for b in x[[bond.a1-1,bond.a2-1]]])
        end
        #
    end every every_
end

#molec(500,10,5)

function plane_js(le)
    syst = load_sdfile(ball_data_path("../test/data/MMFF94-plane.sdf"))
    mm = MMFF94FF(syst)
    x = verlet_at(syst,mm, le, 50, [4])
    at1,at2,at3,at4 = x
    ats(i) = (at1[i],at3[i],at4[i])
    angle(a,b,g) = rad2deg(acos(dot(a-b,g-b)/(norm(a-b)*(norm(g-b)))))
    len = 1000

    for i in 1:le
        trace1 = PlotlyJS.mesh3d(
            x = getproperty.(ats(i),:x),y = getproperty.(ats(i),:y),z = getproperty.(ats(i),:z),
                       colorbar_title="z",showscale=true)
        m = mean([x[1][i], x[3][i], x[4][i]])
        n = cross(x[1][i] - x[3][i], x[1][i] - x[4][i])
        sign_ = -90 <= angle(n,m,x[2][i]-m) <= 90 ? 1 : -1
        no = norm(m-at2[i])
        j = [m, m+sign_*n*no/norm(n)]
        trace2 = PlotlyJS.scatter(x=vcat(getproperty.(j,:x)), y=vcat(getproperty.(j,:y)), z=vcat(getproperty.(j,:z)), type="scatter3d", mode="line")
        trace3 = PlotlyJS.scatter(x = getproperty.(getindex.(x,i),:x),y = getproperty.(getindex.(x,i),:y),z = getproperty.(getindex.(x,i),:z), type="scatter3d", mode="markers")
        p = PlotlyJS.plot([trace1, trace2, trace3])
        PlotlyJS.savefig(p, "temp_gif/im$i.png")
    end
end


function trapeze(a,b,c)
    ab = b-a
    cb = b-c
    trap_x = []
    trap_y = []
    trap_z = []
    #move from c to b in 10% steps
    for v_vert in (c-i*cb for i in 0.0:0.1:1)
        for v_hor in (v_vert - j*ab for j in 0.0:0.1:1)
            push!(trap_x, v_hor.x)
            push!(trap_y, v_hor.y)
            push!(trap_z, v_hor.z)
        end
    end
    return trap_x, trap_y, trap_z
end

#need PlotlyJS
function tors_plot(a, b, c, d)
    f = c - (b-a)
    g = b - (c-d)
    trace1 = PlotlyJS.mesh3d(
        x = getproperty.((b,c,d,g),:x),y = getproperty.((b,c,d,g),:y),z = getproperty.((b,c,d,g),:z),
                   colorbar_title="z",showscale=true)
    trace2 = PlotlyJS.mesh3d(
        x = getproperty.((a,b,c,f),:x),y = getproperty.((a,b,c,f),:y),z = getproperty.((a,b,c,f),:z),
         colorbar_title="z", showscale=true
    )

    angle = rad2deg(acos(dot(a-b,g-b)/(norm(a-b)*(norm(g-b)))))
    anno = attr(
        x=c[1],
        y=c[2],
        z=c[3],
        ax=50,
        ay=0,
        text="Angle=$(trunc(angle))",
        arrowhead=1,
        xanchor="left",
        yanchor="bottom"
    )
    layout = Layout(scene=attr(annotations=[anno]))
    PlotlyJS.plot([trace1, trace2], layout)

end

function look_at_torsions()
    # select a certain timestep in plotly
    # setup: syst =..., mm = ,... x = verlet_at(...)
    ina(idx) = tors_plot(x[1][idx],x[2][idx],x[3][idx],x[4][idx])
    # angles over a verlet_at. For plotting
    syst = load_sdfile(ball_data_path("../test/data/MMFF94-torsion.sdf"))
    mm = MMFF94FF(syst)
    Plots.plot(map((a,b,c,d) -> begin  g = b - (c-d); rad2deg(acos(dot(a-b,g-b)/(norm(a-b)*(norm(g-b))))) end, verlet_at(syst,mm, 2400, 15 ,[2])...))
end


function multi_mole(reps, displace, len, every_, name1, name2; deser=true)

        syst = load_sdfile("$(name1).sdf", Float64)
        syst1 = load_sdfile("$(name2).sdf", Float64)
    

        d = Dict()
        m = Molecule(syst)
        for at in atoms(syst1)
            number  = maximum(atoms_df(syst).number) + 1
            new_at = Atom(m, number, at.element, at.name, at.atom_type, at.r, at.v, at.F,
                at.formal_charge, at.charge, at.radius, at.properties, at.flags)
            d[at.idx] = new_at.idx
        end

        for bond in bonds(syst1)
            Bond(syst, d[bond.a1], d[bond.a2], bond.order, bond.properties, bond.flags)
        end


        for at in eachrow(atoms_df(molecules(syst)[2]))
            at.r += displace
        end

        @show "begin mm"
        mm = MMFF94FF(syst)


    @show "end mm, begin iters"
    

    j = []

    i = 1
    #try
    #    for _ in 1:reps
    #        push!(j, verlet_at(syst,mm,len,10,[]))
    #        println("$(i)/$(reps) done")
    #        i += 1
    #    end
    #finally 
    #    serialize("j_backup_$(i)-$(reps-1)", j)
    #end


    j = [begin println("$(i)/$(reps)"); verlet_at(syst,mm,len,50,[]) end for i in 1:reps]
    h = [vcat((l[i] for l in j)...) for i in 1:natoms(syst)]

    xl = extrema(Iterators.flatten((coords.x for coords in  at) for at in  h))
    xl = (floor(xl[1]), ceil(xl[2]))
    yl = extrema(Iterators.flatten((coords.y for coords in  at) for at in  h))
    yl = (floor(yl[1]), ceil(yl[2]))
    zl = extrema(Iterators.flatten((coords.z for coords in  at) for at in  h))
    zl = (floor(zl[1]), ceil(zl[2]))

    d = Dict(id.idx => i for (i, id) in enumerate(atoms(syst)) )
    @gif for i in 1:length(h[1])
        #scatter([b[i].x for b in x[1]], [b[i].y for b in x[1]], [b[i].z for b in x[1]],
         #   xlim=(-2,2), ylim=(-2,2), zlim=(-2,2))
        Plots.plot(zlim=zl, xlim=xl, ylim=yl, legend=nothing)
        #plot()
        for bond in bonds(syst) 
            Plots.plot!([b[i].x for b in h[[d[bond.a1],d[bond.a2]]]], [b[i].y for b in h[[d[bond.a1],d[bond.a2]]]], [b[i].z for b in h[[d[bond.a1],d[bond.a2]]]])
        end
        #
    end every every_
end

function tra_sim(h,syst, every_,dpi_)
    xl = extrema(Iterators.flatten((coords.x for coords in  at) for at in  h))
    xl = (floor(xl[1]), ceil(xl[2]))
    yl = extrema(Iterators.flatten((coords.y for coords in  at) for at in  h))
    yl = (floor(yl[1]), ceil(yl[2]))
    zl = extrema(Iterators.flatten((coords.z for coords in  at) for at in  h))
    zl = (floor(zl[1]), ceil(zl[2]))

    d = Dict(id.idx => i for (i, id) in enumerate(atoms(syst)) )
    g = @gif for i in 1:length(h[1])
        #scatter([b[i].x for b in x[1]], [b[i].y for b in x[1]], [b[i].z for b in x[1]],
         #   xlim=(-2,2), ylim=(-2,2), zlim=(-2,2))
        Plots.plot(zlim=zl, xlim=xl, ylim=yl, legend=nothing, dpi=dpi_)
        #plot()
        for bond in bonds(syst) 
            Plots.plot!([b[i].x for b in h[[d[bond.a1],d[bond.a2]]]], [b[i].y for b in h[[d[bond.a1],d[bond.a2]]]], [b[i].z for b in h[[d[bond.a1],d[bond.a2]]]])
        end
        #
    end every every_
    return g
end

function sim_multi(syst, reps, len, step, every_;mm=nothing)
    isnothing(mm) && (mm = MMFF94FF(syst);)
    j = []
    try
    #j = [begin println("$(i)/$(reps) en=$(compute_energy(mm))"); verlet_at(syst,mm,len,step,[]) end for i in 1:reps]
        for i in 1:reps    
            push!(j, verlet_at(syst,mm,len,step,[]))
            println("$(i)/$(reps) en=$(compute_energy(mm))")
        end

    catch
        h = [vcat((l[i] for l in j)...) for i in 1:natoms(syst)]
        @show "i died at $(syst.name)"
        return nothing, h

    end
    h = [vcat((l[i] for l in j)...) for i in 1:natoms(syst)]


    #j = [begin println("$(i)/$(reps) en=$(compute_energy(mm))"); verlet_at(syst,mm,len,step,[]) end for i in 1:reps]
    #h = [vcat((l[i] for l in j)...) for i in 1:natoms(syst)]
    @show h[1][1:len:len*reps]

    xl = extrema(Iterators.flatten((coords.x for coords in  at) for at in  h))
    xl = (floor(xl[1]), ceil(xl[2]))
    yl = extrema(Iterators.flatten((coords.y for coords in  at) for at in  h))
    yl = (floor(yl[1]), ceil(yl[2]))
    zl = extrema(Iterators.flatten((coords.z for coords in  at) for at in  h))
    zl = (floor(zl[1]), ceil(zl[2]))

    d = Dict(id.idx => i for (i, id) in enumerate(atoms(syst)) )
    g = nothing
    try
        g = @gif for i in 1:length(h[1])
            #scatter([b[i].x for b in x[1]], [b[i].y for b in x[1]], [b[i].z for b in x[1]],
            #   xlim=(-2,2), ylim=(-2,2), zlim=(-2,2))
            Plots.plot(zlim=zl, xlim=xl, ylim=yl, legend=nothing, dpi=600)
            #plot()
            for bond in bonds(syst) 
                Plots.plot!([b[i].x for b in h[[d[bond.a1],d[bond.a2]]]], [b[i].y for b in h[[d[bond.a1],d[bond.a2]]]], [b[i].z for b in h[[d[bond.a1],d[bond.a2]]]])
            end
            #
        end every every_
        return g,h
    catch 
        @show "i dies at gif"
        return nothing, h
    end
    return g,h
end


function anim_test()
    N=300
    X = LinRange(0, 10, N)
    Y = -3 .+ 7*rand(N)
    Z = +3 .+ 7*rand(N)

    trace = PlotlyJS.scatter(x = [X[1]],  
                    y = [Y[1]],
                    z = [Z[1]],
                    mode="lines",
                    line_width=1.5,
                    line_color="RoyalBlue",
                    type="scatter3d")

    n_frames = length(X)
    frames  = Vector{PlotlyFrame}(undef, n_frames)
    for k in 1:n_frames
        frames[k] = PlotlyJS.frame(data=[attr(x=X[1:k], #update x and y
                                    y=Y[1:k],
                                    z=Z[1:k]
                                    )],
                        layout=attr(title_text="Test frame $k"), #update title
                        name="fr$k", #frame name; it is passed to slider 
                        traces=[0] # this means that the above data update the first trace (here the unique one) 
                            ) 
    end    
  

    ym, yM = extrema(Y)
    zm, zM = extrema(Z)
    layout = Layout(title_text="Test", title_x=0.5,
        width=700, height=450,
                xaxis_range=[-0.1, 10.1], 
                yaxis_range=[ym-1, yM+1],
                zaxis_range=[zm-1, zM+1],
        )
    pl = Plot(trace, layout, frames)
end

function plane_pl(le, displace, del, ev, free)
    syst = load_sdfile(ball_data_path("../test/data/MMFF94-plane2.sdf"))
    at1,at2,at3,at4 = atoms(syst)
    atoms(syst)[free].r += displace
    # m = mean([at1.r, at3.r, at4.r])
    # at1.r += (m-at1.r) * xd
    # at3.r += (m-at3.r) * xd
    # at4.r += (m-at4.r) * xd
    # at2.r += (m-at2.r) * xd2

    mm = MMFF94FF(syst)
    deleteat!(mm.components, del)
    @show length(mm.components)
    x = verlet_at(syst,mm, le, 10, [1,2])
    @show x[1][end]
    @show mean(norm.(x[1].-x[3]))
    @show mean(norm.(x[3].-x[4]))
    @show mean(norm.(x[4].-x[1]))
    at1,at2,at3,at4 = x
    ats(i) = (at1[i],at3[i],at4[i],at1[i])
    angle(a,b,g) = rad2deg(acos(dot(a-b,g-b)/(norm(a-b)*(norm(g-b)))))
    len = 1000
    @show mean([abs(angle(mean([x[1][i], x[3][i], x[4][i]]),x[1][i], x[2][i])) for i in 1:le])
    xl = extrema(Iterators.flatten((coords.x for coords in  at) for at in  x))
    xl = (floor(xl[1]), ceil(xl[2]))
    yl = extrema(Iterators.flatten((coords.y for coords in  at) for at in  x))
    yl = (floor(yl[1]), ceil(yl[2]))
    zl = extrema(Iterators.flatten((coords.z for coords in  at) for at in  x))
    zl = (floor(zl[1]), ceil(zl[2]))


    g = @gif for i in 1:le
        Plots.plot(xlim=xl, ylim=yl, zlim=zl, legend=nothing)
        Plots.plot!([at.x for at in ats(i)], [at.y for at in ats(i)], [at.z for at in ats(i)], linealpha=0.25)
        Plots.scatter!([at.x for at in ats(i)[1:3]],
            [at.y for at in ats(i)[1:3]], [at.z for at in ats(i)[1:3]], msize=6, 
            color="blue")
        Plots.scatter!([x[2][i].x], [x[2][i].y], [x[2][i].z], msize=6, color="red")

        Plots.plot!([x[1][i].x, x[2][i].x],[x[1][i].y, x[2][i].y],[x[1][i].z, x[2][i].z], color="red")
        Plots.plot!([x[3][i].x, x[2][i].x],[x[3][i].y, x[2][i].y],[x[3][i].z, x[2][i].z], color="red")
        Plots.plot!([x[4][i].x, x[2][i].x],[x[4][i].y, x[2][i].y],[x[4][i].z, x[2][i].z], color="red")

        m = mean([x[1][i], x[3][i], x[4][i]])
        n = cross(x[1][i] - x[3][i], x[1][i] - x[4][i])
        sign_ = -90 <= angle(n,m,x[2][i]-m) <= 90 ? 1 : -1
        no = norm(m-at2[i])
        #j = [m, m+sign_*n*no/norm(n)]
        #Plots.plot!([j[1].x, j[2].x], [j[1].y, j[2].y], [j[1].z, j[2].z], color="green", lsize=3)
        #Plots.scatter!([j[1].x, j[2].x], [j[1].y, j[2].y], [j[1].z, j[2].z], color="green",msize=5)

    end every ev

end

function tors_pl(le, displace, del, ev, free)
    syst = load_sdfile(ball_data_path("../test/data/MMFF94-torsion.sdf"))
    at1,at2,at3,at4 = atoms(syst)
    free != 0 && (atoms(syst)[free].r += displace)
    # m = mean([at1.r, at3.r, at4.r])
    # at1.r += (m-at1.r) * xd
    # at3.r += (m-at3.r) * xd
    # at4.r += (m-at4.r) * xd
    # at2.r += (m-at2.r) * xd2

    mm = MMFF94FF(syst)
    deleteat!(mm.components, del)
    @show length(mm.components)
    x = verlet_at_stat(syst,mm, le, 10, [1,2], free)
    @show x[1][end]
    at1,at2,at3,at4 = x
    ats(i) = (at1[i],at2[i],at3[i],at4[i])
    angle(a,b,g) = rad2deg(acos(dot(a-b,g-b)/(norm(a-b)*(norm(g-b)))))
    len = 1000
    angles = [abs(angle(mean([x[1][i], x[3][i], x[4][i]]),x[1][i], x[2][i])) for i in 1:le]
    xl = extrema(Iterators.flatten((coords.x for coords in  at) for at in  x))
    xl = (floor(xl[1]), ceil(xl[2]))
    yl = extrema(Iterators.flatten((coords.y for coords in  at) for at in  x))
    yl = (floor(yl[1]), ceil(yl[2]))
    zl = extrema(Iterators.flatten((coords.z for coords in  at) for at in  x))
    zl = (floor(zl[1]), ceil(zl[2]))


    @gif for i in 1:le
        Plots.plot(xlim=xl, ylim=yl, zlim=zl, legend=nothing)
        Plots.plot!([at.x for at in ats(i)], [at.y for at in ats(i)], [at.z for at in ats(i)], linealpha=1)
        Plots.scatter!([at.x for at in ats(i)], [at.y for at in ats(i)], [at.z for at in ats(i)], linealpha=1, color="blue", msize=6)
    end every ev
    return x
end

function sep(names)
    name = first(names)
    load_op = endswith(name,".sdf") ? load_sdfile : load_pdb
    sys = load_op(name)
    d = Dict()

    for name in names[2:end]
        load_op = endswith(name,".sdf") ? load_sdfile : load_pdb
        new_sys = load_op(name)
        m = Molecule(sys)

        for at in atoms(new_sys)
            number  = maximum(atoms_df(sys).number) + 1
            new_at = Atom(m, number, at.element, at.name, at.atom_type, at.r, at.v, at.F,
                at.formal_charge, at.charge, at.radius, at.properties, at.flags)
            d[at.idx] = new_at.idx
        end

        for bond in bonds(new_sys)
            Bond(sys, d[bond.a1], d[bond.a2], bond.order, bond.properties, bond.flags)
        end

    end 
    return sys
end

function central(sys::System{T}, dist) where T <:Real
    cent = molecules(sys)[1]
    mean_c = Vector3{T}(mean(atoms_df(cent).r))
    rotz(a) = [cos(deg2rad(a)) -sin(deg2rad(a)) 0; sin(deg2rad(a)) cos(deg2rad(a)) 0; 0 0 1]
    mov(a,b) = b-a

    rot = 360/(length(molecules(sys))-1)
    for (i,m) in enumerate(molecules(sys)[2:end])
        mean_m = Vector3{T}(mean(atoms_df(m).r))
        #move molecule onto central molecule
        atoms_df(m).r = map(at-> at + mov(mean_m, mean_c), atoms_df(m).r)
        #move molecule a certain dist away
        atoms_df(m).r = map(at-> at + Vector3{T}(dist, 0, 0), atoms_df(m).r)
        #move molecule to origin
        atoms_df(m).r = map(at-> at + mov(mean_c, Vector3{T}(0,0,0)), atoms_df(m).r)
        #rotate
        atoms_df(m).r = map(at -> Vector3{T}(rotz(i*rot) *at) , atoms_df(m).r)
        #move molec back around central mol
        atoms_df(m).r = map(at-> at + mean_c, atoms_df(m).r)
    end

end

function glucose_water(;cops,dist,reps,len,step,every_)
    sys = sep(vcat("water.sdf",repeat(["water.sdf"],cops)))
    central(sys,dist)#
    mm = MMFF94FF(sys)
    return sim_multi(sys, reps, len, step, every_;mm=mm)

end

function test_()
    names = ["1p8b","3i40","4ey0"]
    for name in names
        @show "doing $name"
        sys = load_sdfile("$(name).sdf")
        mm = MMFF94FF(sys)
        gif, h = sim_multi(mm.system, 1000, 50, 35, 100,mm=mm)
        serialize("$(name)_mm.ser", mm)
        serialize("$(name)_tra.ser", h)

    end

end

function water_dipole(dist)
    sys = load_sdfile("glucose.sdf")
    new_sys = load_sdfile("water.sdf")
    m = Molecule(sys)
    d = Dict()

    for at in atoms(new_sys)
        number  = maximum(atoms_df(sys).number) + 1
        new_at = Atom(m, number, at.element, at.name, at.atom_type, at.r, at.v, at.F,
            at.formal_charge, at.charge, at.radius, at.properties, at.flags)
        d[at.idx] = new_at.idx
    end

    for bond in bonds(new_sys)
        Bond(sys, d[bond.a1], d[bond.a2], bond.order, bond.properties, bond.flags)
    end

    for at in atoms(m)
        at.r += dist
    end
    mm = MMFF94FF(sys)
    deleteat!(mm.components, [2])
    return mm
    
end