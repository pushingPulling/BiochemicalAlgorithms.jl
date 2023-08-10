
@auto_hash_equals mutable struct MTorsion{T<:Real}
    type::Int
    at1::Atom{T}
    at2::Atom{T}
    at3::Atom{T}
    at4::Atom{T}
    v1::T
    v2::T
    v3::T
    energy::T
    angle::T
    heuristic::Bool
    
end

mutable struct MTorsionComponent{T<:Real} <: AbstractForceFieldComponent{T}
    ff::ForceField{T}
    torsions::Vector{MTorsion{T}}
    energy::Dict{String, T}
    rings_cache::Tuple{Vector{Vector{Atom}}, Vector{Vector{Atom}}}
end



function MTorsionComponent{T}(ff, all_rings, aromatic_rings) where {T<:Real}
    all_rings_ = collect(Iterators.flatten(all_rings))
    aromatic_rings_ = collect(Iterators.flatten(aromatic_rings))
    comp = MTorsionComponent{T}(ff, MTorsion{T}[], Dict{String, T}(),
        (all_rings_, aromatic_rings_))

    setup!(comp)
    update!(comp)
    comp
end

function setup!(mtc::MTorsionComponent{T}) where {T<:Real}
    equivs = mtc.ff.parameters.sections["Equivalences"].data
    torsion_df = mtc.ff.parameters.sections["Torsions"].data
    at_type_props = mtc.ff.parameters.sections["AtomTypeProperties"].data
    all_rings = mtc.rings_cache[1]
    aro_rings = mtc.rings_cache[2]
    gdf_torsion = groupby(torsion_df, [:at, :I, :J, :K, :L])
    gdf_l = groupby(non_hydrogen_bonds_df(mtc.ff.system), [:a1])
    gdf_r = groupby(non_hydrogen_bonds_df(mtc.ff.system), [:a2])

    empty!(mtc.torsions)
    empty!(mtc.energy)

    #gets all bonds in which at is a part of; the grouped dataframes version
    _get_bonds(at) = eachrow(vcat(filter(!isnothing, (get(gdf_l, (a1=at.idx,), nothing),
        get(gdf_r, (a2=at.idx,), nothing)))...))

    #starts at line 141
    #each  at_x = at
    sys = mtc.ff.system

    #find possible torisons with this pattern (at = atom, nhb = nonhydrogenbond)
    # |at1| <--nhb2--> |at2| <--nhb1--> |at3| <--nhb3--> |at4|
    for at in atoms(sys)        #this atom is |at2| in the pattern above
        #non hydrobond at_x->h1
        for nhb1 in _get_bonds(at)
            at.idx != nhb1.a1 && continue       #ensures we use a bond only once
            at2 = at
            at3_orig = atom_by_idx(sys, nhb1.a2)#get_partner(nhb1, at2)
            #need to copy this value here, because it can get swapped down there. Need the original at3 at the start of this loop
            for nhb2 in _get_bonds(at)#[cur_bond + 1 : end] #non_hydrobond at_x->h2
                nhb1.idx == nhb2.idx && continue
                at1_orig = get_partner(nhb2, at)
                for nhb3 in _get_bonds(at3_orig)
                    #restore values before swap
                    at3 = at3_orig# deepcopy(at3_orig)
                    at2 = at#deepcopy(at)
                    at4 = get_partner(nhb3, at3)
                    at1 = at1_orig#deepcopy(at1_orig)


                    (at4.idx == at2.idx || at4.idx == at1.idx) && continue
                    
                    at2.atom_type != "0" && at_type_props[at_idx(at2.atom_type), :lin] |> Bool && continue
                    at3.atom_type != "0" && at_type_props[at_idx(at3.atom_type), :lin] |> Bool && continue

                    #original BALL code checks if there is a ring of size 3 with all atoms in it
                    # that does never happen, as at1 != at2 != at3 != at4 and a ring of size 3 can't hold 4 atoms

                    if (parse(Int8, at2.atom_type) > parse(Int8, at3.atom_type)) ||
                        (parse(Int8, at2.atom_type) == parse(Int8, at3.atom_type) && parse(Int8, at1.atom_type) > parse(Int8, at4.atom_type))
                        at1, at2, at3, at4 = at4, at3, at2, at1
                    end


                    torsion_type = get_torsion_type(at1, at2, at3, at4,
                        nhb2, nhb1, nhb3, all_rings, aro_rings)


                    torsion_v1, torsion_v2, torsion_v3 = get_torsion_params(
                        torsion_type, at1, at2, at3, at4, gdf_torsion, equivs
                    )
                    heuristic_used = false
                    if ismissing(torsion_v1)
                        torsion_v1, torsion_v2, torsion_v3 = calc_torsion_heuristic(
                            at2, at3, nhb1, at_type_props, aro_rings
                        )
                        heuristic_used = true

                        if ismissing(torsion_v1)
                            push!(mtc.ff.unassigned_atoms, at1, at2, at3, at4)


                            @warn """Could not find Torsion params for atoms 
                                 $(at1.name):$(at1.atom_type), 
                                 $(at2.name):$(at2.atom_type), 
                                 $(at3.name):$(at3.atom_type), 
                                 $(at4.name):$(at4.atom_type), 
                            """
                        end
                    end
                    all(ismissing, (torsion_v1, torsion_v2, torsion_v3)) && continue
                    push!(mtc.torsions, MTorsion{T}(torsion_type, at1, at2,
                        at3, at4, torsion_v1, torsion_v2, torsion_v3, zero(T),
                        zero(T), heuristic_used)
                    )
                end
            end
        end
    end
end

function update!(mtc::MTorsionComponent{T}) where {T<:Real}
    nothing
end

function compute_energy(mtc::MTorsionComponent{T}) where {T<:Real}
    mapreduce(compute_energy, +, mtc.torsions, init = zero(T))
end

function compute_energy(torsion::MTorsion{T}) where {T<:Real}
    a21 = torsion.at1.r - torsion.at2.r
    a23 = torsion.at3.r - torsion.at2.r
    a34 = torsion.at4.r - torsion.at3.r

    cross2321 = normalize(cross(a23, a21))
    cross2334 = normalize(cross(a23, a34))

    any(isnan, first.((cross2321, cross2334))) && return 0

    cosphi = dot(cross2321, cross2334)
    cosphi = max(min(cosphi, 1), -1)
    phi = acos(cosphi)

    e = k_torsion * (torsion.v1 * (1 + cosphi) 
        + torsion.v2 * (1 - cos(2 * phi))
        + torsion.v3 * (1 + cos(3 * phi)))

end

function compute_forces(torsion::MTorsion{T}) where {T<:Real}
    angle(a::Vector3{T}, b::Vector3{T}) = acos(max(min(dot(a,b) / (sqrt(norm(a) * norm(b))),1),-1))

    at1 = torsion.at1
    at2 = torsion.at2
    at3 = torsion.at3
    at4 = torsion.at4

    ij = at2.r - at1.r
    jk = at3.r - at2.r
    kl = at4.r - at3.r
    ij_l = norm(ij)
    jk_l = norm(jk)
    kl_l = norm(kl)

    any(==(0), (ij_l, jk_l, kl_l)) && return 0
    
    ijk_a = angle(ij, jk)
    jkl_a = angle(jk, kl)

    ij /= ij_l
    jk /= jk_l
    kl /= kl_l

    sin_j = sin(ijk_a)
    sin_k = sin(jkl_a)

    rsj = inv(ij_l * sin_j^2)
    rsk = inv(kl_l * sin_k^2)

    rrcj = (ij_l / jk_l) * (-cos(ijk_a))
    rrck = (kl_l / jk_l) * (-cos(jkl_a))

    a = cross(ij, jk)
    b = cross(jk, kl)
    c = cross(a, b)

    d1 = dot(c, jk)
    d2 = dot(a, b)
    angle_d1d2 = atan(d1, d2)

    di = -a * rsj
    dl =  b * rsk
    dj = di * (rrcj -1) - dl * rrck
    dk = -(di + dj + dl)

    c1 = (torsion.v1 * sin(angle_d1d2) - 2 * torsion.v2 * sin(angle_d1d2 * 2) +
        3 * torsion.v3 * sin(angle_d1d2 * 3) )

    c1 *= 0.5u"cal" |> u"J" |> ustrip
    c1 *= 1u"kJ/mol/Å"/Constants.N_A |> u"N" |> ustrip

    torsion.at1.F += di * c1
    torsion.at2.F += dj * c1
    torsion.at3.F += dk * c1
    torsion.at4.F += dl * c1 
    

end

function compute_forces(mtc::MTorsionComponent{T}) where {T<:Real}
    map(compute_forces, mtc.torsions)
    nothing
end

function calc_torsion_heuristic(at_j, at_k, bond, at_props, aro_rings)

    v1 = 0
    v2 = 0
    v3 = 0
    (at_j.atom_type == "0" || at_k.atom_type == "0") && return missing, missing, missing
    props_j = at_props[at_idx(at_j.atom_type), :]
    props_k = at_props[at_idx(at_j.atom_type), :]

    #rule a) Linearity
    
    (props_j.lin == 1 || props_k.lin == 1) && return v1, v2, v3


    uj = _get_U(at_j)
    uk = _get_U(at_k)
    vj = _get_V(at_j)
    vk = _get_V(at_k)

    beta = 6
    n = (props_j.crd - 1) * (props_k.crd - 1)



    #rule b) both atomm types are aromatic
    if any(r -> at_k in r && at_j in r, aro_rings)
        (uj == -1 || uk == -1) && return missing, missing, missing

        #one trivalent and one tetravalent?
        beta = props_j.val * props_k.val == 12 ? (beta = 3) : (beta = 6)
        beta = props_j.pilp * props_k.pilp == 12 ? (l = 0.3) : (l = 0.5)

        v2 = beta * l * sqrt(uj * uk)
        return v1, v2, v3
    end

    #rule c) double bond


    if bond.order == BondOrder.Double
        (uj == -1 || uk == -1) && return missing, missing, missing
        props_j.mltb * props_k.mltb == 4 ? (l = 1) : (l = 0.4)

        v2 = beta * l * sqrt(uj * uk)
        return v1, v2, v3
    end


    #rule d) both atoms must be trtracoordinate
    if props_j.crd * props_k.crd == 16
        (vj == -1 || vk == -1) && return missing, missing, missing
        v3 = sqrt(vj * vk) / n

        return v1, v2, v3
    end

    #rule e) and f) j or k tetravalent
    props = (props_j, props_k)
    for (i, at_type) in enumerate(props)
        if props[mod1(i + 1, 2)].crd == 4     # use the opposite props in this line.
            at_type.crd == 3 && (at_type.val == 4 || at_type.val == 34 || at_type.mltb) &&
                return v1, v2, v3

            at_type.crd == 2 && (at_type.val == 3 || at_type.mltb) &&
                return v1, v2, v3

            (vj == -1 || vk == -1) && return missing, missing, missing
            v3 = sqrt(vj * vk) / n

            return v1, v2, v3
        end
    end

    #rule g) k trivalent
    if (bond.order == BondOrder.Single && props_j.mltb == 1 && props_k.mltb == 1) ||
       (props_j.pilp == 1 && props_k.mltb == 1) ||
       (props_j.pilp == 1 && props_k.mltb == 1)

        props_j.pilp == 1 && props_k.pilp == 1 && return v1, v2, v3
        (uk == -1 || uj == -1) && return missing, missing, missing

        if (props_j.pilp == 1 && props_k.mltb == 1)
            if props_j.mltb == 1
                l = 0.5
            elseif at_j.element == Elements.H && at_k.element == Elements.H
                l = 0.3
            end
        elseif (props_k.pilp == 1 && props_j.mltb == 1)
            if (props_k.mltb == 1)
                l = 0.5
            elseif at_j.element == Elements.H && at_k.element == Elements.H
                l = 0.3
            end
        elseif (props_j.mltb == 1 || props_k.mltb == 1) &&
               (at_j.element != Elements.C || at_k.element != Elements.C)

            l = 0.4
        else
            l = 0.15        #default value
        end

        v2 = beta * l * sqrt(uj * uk)
        return v1, v2, v3

    end

    #rule h) saturated centers, at most trivalent

    if (at_j.element == Elements.O || at_j.element == Elements.S) &&
       (at_k.element == Elements.O || at_k.element == Elements.S)

        es = Int8(at_j.element) + Int8(at_k.element)

        if es == 24
            v2 = -4
        elseif es == 32
            v2 = -8
        else
            v2 = -2
        end

        return v1, v2, v3

    end

    #euqation 22
    (vk == -1 || vj == -1) && return missing, missing, missing

    v3 = sqrt(vj * vk) / n
    return v1, v2, v3
end


function _get_U(at)
    idx = Int8(at.element)
    #values taken from paper IV
    5 < idx < 9 && return 2.0
    13 < idx < 17 && return 1.25
    #values taken from CHARMM implementation
    31 < idx < 35 && return 2.0
    49 < idx < 53 && return 2.0
    idx < 86 && return 0.1

    return -1

end

function _get_V(at)
    idx = Int8(at.element)
    # values taken from paper IV
    idx == 6 && return 2.12
    idx == 7 && return 1.50
    idx == 8 && return 0.20
    idx == 14 && return 1.22
    idx == 15 && return 2.40
    idx == 16 && return 0.49
    # values tk&& en from the CHARMM implementation:
    idx == 32 && return 0.7
    idx == 33 && return 1.5
    idx == 34 && return 0.34
    idx == 50 && return 0.2
    idx == 51 && return 1.1
    idx == 52 && return 0.3
    idx < 86 && return 0.3

    return -1
end


function get_torsion_params(torsion_type, at1, at2, at3, at4, gdf, equivs)
    at1_t, at2_t, at3_t, at4_t = getproperty.((at1, at2, at3, at4), :atom_type)

    row = get(gdf, (at=torsion_type, I=              at1_t      , J=at2_t, K=at3_t, L=              at4_t), missing)
    !ismissing(row) && return row[:V1, :V2, :V3]
    any(==("0"), (at1_t, at4_t)) && return missing, missing, missing

    row = coalesce(
        
        get(gdf, (at=torsion_type, I=equivs[at_idx(at1_t), :b], J=at2_t, K=at3_t, L=equivs[at_idx(at4_t), :b]), missing),
        get(gdf, (at=torsion_type, I=equivs[at_idx(at1_t), :c], J=at2_t, K=at3_t, L=equivs[at_idx(at4_t), :c]), missing),
        get(gdf, (at=torsion_type, I=equivs[at_idx(at1_t), :d], J=at2_t, K=at3_t, L=equivs[at_idx(at4_t), :d]), missing),
        get(gdf, (at=torsion_type, I=equivs[at_idx(at1_t), :e], J=at2_t, K=at3_t, L=equivs[at_idx(at4_t), :e]), missing)
    )

    !ismissing(row) && (return row[:V1, :V2, :V3])
    return missing, missing, missing

end

function get_torsion_type(at1, at2, at3, at4, bond_1_2, bond_2_3, bond_3_4,
         all_rings, aro_rings)

    # if these atoms are in a ring of size 4 or in 2 rings of size 3 return 4
    any(r -> length(r) == 4 && all(at -> at in r, (at1, at2, at3, at4)),
            (r for r in all_rings if length(r) == 4)
        ) &&
        at3 in bonds(at1) &&
        at4 in bonds(at2) &&
        return 4

    if bond_2_3.order == BondOrder.Single &&
       !any(r -> all(at -> at in r, (at2, at3)), aro_rings)

        hasproperty(bond_2_3, :MMFF94SBMB) && return 1

        (
            hasproperty(bond_1_2, :MMFF94SBMB) ||
            hasproperty(bond_3_4, :MMFF94SBMB)
        ) && return 2
    end


    any(r -> length(r) == 5 && all(at -> at in r, (at1, at2, at3, at4)),
            (r for r in all_rings if length(r) == 4)
        ) &&
        any(at -> at.atom_type == "1", (at1, at2, at3, at4)) &&
        return 5

    #default is 0
    return 0

end