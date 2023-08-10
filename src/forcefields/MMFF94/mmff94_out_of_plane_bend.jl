@auto_hash_equals mutable struct MOutOfPlaneBend{T<:Real}
    i::Atom{T}
    j::Atom{T}
    k::Atom{T}
    l::Atom{T}
    koop::T
end

@auto_hash_equals mutable struct MOutOfPlaneComponent{T<:Real} <: AbstractForceFieldComponent{T}

    name::String
    ff::ForceField{T} 
    cache::Dict{Symbol, Any}
    energy::Dict{String, T}
    out_of_plane_bends::Vector{MOutOfPlaneBend}


    function MOutOfPlaneComponent{T}(ff::ForceField{T}) where {T<:Real}
        this = new("moop", ff, Dict{Symbol, Any}(), Dict{String, T}(),
            MOutOfPlaneBend{T}[]
        )

        setup!(this)
        update!(this)

        this
    end
end

function update!(moop::MOutOfPlaneComponent{T}) where T<:Real
    nothing
end

function setup!(moop::MOutOfPlaneComponent{T}) where T<:Real
    # check if enabled bla bla
    empty!(moop.cache)
    empty!(moop.energy)
    empty!(moop.out_of_plane_bends)
    # naive method of iterating through all atoms and picking out the bonds
    # is 10x slower than this approach (on a system with 150 atoms)
    df = non_hydrogen_bonds_df(moop.ff.system)
    isempty(df) && return
    gdf_l = groupby(df, [:a1])
    gdf_r = groupby(df, [:a2])
    gdf_oop = groupby(moop.ff.parameters.sections["OutOfPlane"].data , [:a, :b, :c, :d])
    equivs = moop.ff.parameters.sections["Equivalences"].data
    at_idxs = (minimum(minimum.((df.a1, df.a2))), maximum(maximum.((df.a1,df.a2))))
    #makes it possible to concatenate an empty DataFrame with a populated one
    dummy_df = DataFrame(propertynames(df) .=> (Int64[], Int64[], Int64[], BondOrder.T[], Dict{Symbol, Any}[], Set{Symbol}[]))

    for index in at_idxs[1]:at_idxs[2]
        f = vcat(get(gdf_l, (a1=index,), dummy_df), get(gdf_r, (a2=index,),dummy_df))
        size(f)[1] != 3 && continue;
        #get the atoms in f
        ats = atom_by_idx.((moop.ff.system for _ in size(f)[1]), (a.a1 == index ? a.a2 : a.a1 for a in eachrow(f)))
        
        central_at_idx = index
        central_at = atom_by_idx(moop.ff.system, index)
        j = central_at
        i, k, l = sort(ats, by=at -> parse(Int8, at.atom_type))

        koop = get_oop_params(i.atom_type, j.atom_type, k.atom_type, l.atom_type,
            gdf_oop, equivs)

        if ismissing(koop)
            @warn """Cannot find find OOP params for at types $(i.atom_type),\
            $(j.atom_type), $(k.atom_type), $(l.atom_type)"""
            #assign into unassigned atoms
            continue
        end

        koop != 0.0 && push!(moop.out_of_plane_bends, MOutOfPlaneBend{T}(i, j, k, l, koop))
    end
end

function get_oop_params(a1_t, a2_t, a3_t, a4_t, gdf, equivs)

    row = get(gdf, (a=       at_idx(a1_t)     , b=at_idx(a2_t), c=       at_idx(a3_t) , d=       at_idx(a4_t)  ), missing)
    !ismissing(row) && return only(row.koop)
    any(==("0"), (a1_t, a3_t, a4_t)) && return missing
    row = coalesce(
        
        get(gdf, (a=equivs[at_idx(a1_t), :b], b=at_idx(a2_t), c=equivs[at_idx(a3_t), :b], d=equivs[at_idx(a4_t), :b] ), missing),
        get(gdf, (a=equivs[at_idx(a1_t), :c], b=at_idx(a2_t), c=equivs[at_idx(a3_t), :c], d=equivs[at_idx(a4_t), :c] ), missing),
        get(gdf, (a=equivs[at_idx(a1_t), :d], b=at_idx(a2_t), c=equivs[at_idx(a3_t), :d], d=equivs[at_idx(a4_t), :d] ), missing),
        get(gdf, (a=equivs[at_idx(a1_t), :e], b=at_idx(a2_t), c=equivs[at_idx(a3_t), :e], d=equivs[at_idx(a4_t), :e] ), missing),
    )

    !ismissing(row) && return only(row.koop)
    
    return missing
    
end

function compute_energy(moop::MOutOfPlaneComponent{T}) where T<:Real
    mapreduce(compute_energy, + , moop.out_of_plane_bends, init = zero(T))
end

function compute_forces(moop::MOutOfPlaneComponent{T}) where T<:Real
    map(compute_forces, moop.out_of_plane_bends)
    nothing
end

function compute_energy(bend::MOutOfPlaneBend{T}) where T<:Real
    vs_0 = normalize(Vector3{Float64}(bend.j.r - bend.i.r))
    vs_1 = normalize(Vector3{Float64}(bend.j.r - bend.k.r))
    vs_2 = normalize(Vector3{Float64}(bend.j.r - bend.l.r))

    n_0 =  normalize(cross(vs_1, vs_2))
    n_1 =  normalize(cross(vs_0, vs_2))
    n_2 =  normalize(cross(vs_0, vs_1))

    any(isnan, first.((vs_0, vs_1, vs_2, n_0, n_1, n_2))) && return 0

    k = Float64(0.021922) * bend.koop        #k0 constant = 0.043844/2 = 0.021922
    e = Float64(0.)

    for (v, n) in ((vs_0, n_0), (vs_1, n_1), (vs_2, n_2))
        #intersection_angle::Float64 = rad2deg(asin(abs(dot(v, n))))
        intersection_angle::Float64 = dot(v,n) |> abs |> asin |> rad2deg
        e += k * intersection_angle^2
    end

    e = e * joule_per_cal
    return e

end

function compute_forces(bend::MOutOfPlaneBend{T}) where T<:Real
    FC = 0.043844 * (180/Constants.pi)^2

    center_at = bend.j
    partners = [bend.i, bend.k, bend.l]

    deltas = SVector{3, Float64}[bend.i.r - center_at.r,
        bend.k.r - center_at.r, bend.l.r - center_at.r]
    lens = norm.(deltas)
    deltas ./= lens
    any(==(0), lens) && return 0

    ##
    for permutation in ([1,2,3], [1,3,2], [2,3,1])
        i,      k,      l       = partners[permutation]
        ji,     jk,     jl      = deltas[permutation]
        len_ji, len_jk, len_jl  = lens[permutation]

        an = cross(ji, jk)
        bn = cross(jk, jl)
        cn = cross(jl, ji)

        cos_theta = dot(ji, jk)
        theta = acos(cos_theta)
        theta == 0 && continue

        sin_dl = dot(an, jl) / sin(theta)
        dl = asin(sin_dl)
        dl == 0 && continue

        cos(dl) < 0.0001 && continue     #if wilson angle is equal to 9 degree

        c1 = -dl * FC * bend.koop * force_prefactor * joule_per_cal
        tmp = cos(dl) / c1

        d_l = ((an / sin(theta) - jl * sin(dl)) / len_jl) / tmp
        d_i = (((bn + (((-ji + jk * cos(theta)) * sin(dl)) / sin(theta))) / len_ji) / tmp) / sin(theta)
        d_k = (((cn + (((-jk + ji * cos(theta)) * sin(dl)) / sin(theta))) / len_jk) / tmp) / sin(theta)

        i.F += d_i
        k.F += d_k
        l.F += d_l
        center_at.F -= (d_i + d_k + d_l)

    end

end