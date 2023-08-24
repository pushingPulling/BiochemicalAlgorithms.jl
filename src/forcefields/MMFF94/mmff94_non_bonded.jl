mutable struct ESConstantTerms{T<:Real}
    ac::T
    bc::T
    c3::T
    d5::T
    denom::T
    add::T
    addr::T
    const_::T
    constr_::T
    CONSTANT_::T #in Joule

    function ESConstantTerms{T}(cn::T, cf::T) where {T<:Real}    #cn = cut on; cf = cut off
        cn += 0.05
        cf += 0.05
        denom = (cf^2-cn^2)^3
        ac = cf^4* (cf^2-3*cn^2) / denom
        bc = (6*cn^2 * cf^2) / denom
        c3 = -(cn^2 + cf^2) / denom
        d5 = (2/denom) / 5
        const_ = bc * cf - ac/cf + c3*cf^3 + d5 * cf^5
        add = (cf^2 * cn^2 * (cf - cn) - (cf^5 - cn^5)/5) * 8/denom
        constr_ = 2 * bc *log(cf) - ac/cf^2 + 3 * c3 * cf^2 + cf^4/denom
        CONSTANT_ = 1389.3875744
        addr = (12 * cn^2 * cf^2 * log(cf/cn) - 3 * (cf^4 - cn^4)) / denom

        new(ac, bc, c3, d5, denom, add, addr, const_, constr_, CONSTANT_)
    end
    
end

# The buffered potential takes two Parameters (for our purpose M=14, N=7)
@auto_hash_equals mutable struct BufferedVdWInteraction{T<:Real,M,N}
    distance::T
    a1::Atom{T}
    a2::Atom{T}
    rij::T
    rij_7::T
    eij::T
end

# ES interactions depend on a certain cutoff/cuton C. 
# This is going to simplify calculations and hopefully speed up the code
# ddd is distance dependelt dielectric
@auto_hash_equals mutable struct BufferedESInteraction{T<:Real, M, N}
    distance::T
    cut_off::T
    cut_on::T
    a1::Atom{T}
    a2::Atom{T}
    qi::T
    qj::T
    is_1_4::Bool
    ddd::Bool
    es_switch::Bool
    es_constants::ESConstantTerms{T}
end

@auto_hash_equals mutable struct MNonBondedComponent{T<:Real} <: AbstractForceFieldComponent{T}

    name::String
    ff::ForceField{T} 
    cache::Dict{Symbol, Any}
    energy::Dict{String, T}
    vdw_interactions::AbstractVector{BufferedVdWInteraction{T, 14, 7}}
    es_interactions::AbstractVector{BufferedESInteraction{T, 14, 7}}
    vdw_cache::Dict{Vector{String}, NTuple{3,T}}
    es_switch::Bool
    es_constants::ESConstantTerms{T}

    function MNonBondedComponent{T}(ff::ForceField{T}) where {T<:Real}
        this = new("mnb", ff, Dict{Symbol, Any}(), Dict{String, T}(),
            BufferedVdWInteraction{T, 14, 7}[], BufferedESInteraction{T, 14, 7}[],
            Dict{Vector{String}, NTuple{3,T}}(), false,
            ESConstantTerms{T}(ff.options[:electrostatic_cuton],
                ff.options[:electrostatic_cutoff])
        )

        setup!(this)
        update!(this)

        this
    end
end

#a struct for values that only depend on the ES_cutoff and ES_cuton


function get_VdW_params(at1_t, at2_t, vdw_df, mnb::MNonBondedComponent{T}) where {T<:Real}

    cache = get(mnb.vdw_cache, sort([at1_t, at2_t]), missing)
    !ismissing(cache) && return cache

    (at1_t == "0" || at2_t == "0") && return missing, missing, missing
    row1 = vdw_df[at_idx(at1_t), :]
    row2 = vdw_df[at_idx(at2_t), :]

    ri = row1.a_i * row1.al_i^(1//4)
    rj = row2.a_i * row2.al_i^(1//4)

    l = (ri-rj) / (ri+rj)

    if at1_t != at2_t
        rij = (row1.da == "D" || row2.da =="D") ? 
            0.5 * (ri + rj) :
            0.5 * (ri + rj) * (1. + 0.2 * (1. - exp(-12. * l^2 )))
    else
        rij = ri
    end

    up = 181.16 * row1.g_i * row2.g_i * row1.al_i * row2.al_i
    lo = sqrt(row1.al_i / row1.n_i) + sqrt(row2.al_i / row2.n_i)

    eij = (up/lo) / rij^6

    if at1_t != at2_t && (row1.da, row2.da) in (("D", "A"), ("A", "D"))
        rij *= 0.8
        eij *= 0.5
    end

    mnb.vdw_cache[sort([at1_t, at2_t])] = (rij, rij^7, eij)

end

function setup!(mnb::MNonBondedComponent{T}) where T<:Real
    #i am assuming es and vdw are enabled 
    #these properties dont actually exist in the options field
    #else: !haskey(mnb.ff.options[:es_enabled]) && haskey(mnb.ff.options[:vdw_enabled]) && continue

    #setting values from options
    nonbonded_cutoff = get(mnb.ff.options, :nonbonded_cutoff, 0.)
    vdw_cut_off      = get(mnb.ff.options, :vdw_cutoff, 0.)
    vdw_cut_on       = get(mnb.ff.options, :vdw_cuton, 0.)
    es_cut_off       = get(mnb.ff.options, :electrostatic_cutoff, 0.)
    es_cut_on        = get(mnb.ff.options, :electrostatic_cuton, 0.)
    ddd              = get(mnb.ff.options, :distance_dependent_dielectric, false)

    periodic_box = [
        mnb.ff.options[:periodic_box_width], 
        mnb.ff.options[:periodic_box_height],
        mnb.ff.options[:periodic_box_depth]
    ]

    if get(mnb.ff.options, :periodic_boundary_conditions, false)
        max_cut_off = 0.5 * min(periodic_box)
        if cur_off > max_cut_off
            @error """MMFF94NonBonded::setup(): Electrostatic cutoff may not 
            exceed half the box dimension when using periodic boundary conditions! 
            Aborting setup."""
            return
        end
    end

    #271
    if es_cut_off > 0
        if es_cut_off > es_cut_on && nonbonded_cutoff >= es_cut_off
            mnb.es_switch = true    #this ties to the es_constants struct
        else
            mnb.es_switch = false
            @error "MMFF94NonBonded::setup(): Invalid ES cutof/cutoff values!"
        end  
    end

    if vdw_cut_off > 0
        if !(vdw_cut_off > vdw_cut_on && nonbonded_cutoff >= vdw_cut_off)
            @error "MMFF94NonBonded::setup(): Invalid VdW cutof/cutoff values! Aborting setup."
        end
    end

    #taken from AMBER implementation. vicinal relationship is needed later
    mg = MetaGraph(non_hydrogen_bonds_df(mnb.ff.system), :a1, :a2)

    nh_0 = Set.(map(v -> mg.vprops[v][:name], vertices(mg)))

    to_set(nh) = Set{Pair{Int64, Int64}}(
        a => b for i in eachindex(nh_0) for a in nh_0[i] for b in nh[i]
    )

    compute_neighborhood(level) = 
        map(v->Set(map(v->mg.vprops[v][:name], neighborhood(mg, v, level))), vertices(mg))

    nh_1 = compute_neighborhood(1)
    nh_2 = compute_neighborhood(2)
    nh_3 = compute_neighborhood(3)

    bond_cache    = to_set(setdiff.(nh_1, nh_0))
    geminal_cache = to_set(setdiff.(nh_2, nh_1))
    vicinal_cache = to_set(setdiff.(nh_3, nh_2))

    # remember those parts that stay constant when only the system is updated
    mnb.cache[:nonbonded_cutoff]                = nonbonded_cutoff
    mnb.cache[:periodic_box]                    = periodic_box
    mnb.cache[:electrostatic_cutoff]            = es_cut_off
    mnb.cache[:electrostatic_cuton]             = es_cut_on
    mnb.cache[:vicinal_cache]                   = vicinal_cache
    mnb.cache[:bond_cache]                      = bond_cache
    mnb.cache[:geminal_cache]                   = geminal_cache
    mnb.cache[:distance_dependent_dielectric]   = ddd

    nothing

end

function update!(mnb::MNonBondedComponent{T}) where T<:Real


    # calc nonbonded atom pairs.
    # change some params using vdw database
    periodic_box     = mnb.cache[:periodic_box]
    nonbonded_cutoff = mnb.cache[:nonbonded_cutoff]
    es_cut_off       = mnb.cache[:electrostatic_cutoff]
    es_cut_on        = mnb.cache[:electrostatic_cuton]
    ddd              = get(mnb.ff.options, :distance_dependent_dielectric, false)
    vdw_df           = mnb.ff.parameters.sections["VanDerWaals"].data

    neighbors = ((mnb.ff.options[:periodic_boundary_conditions]) 
        ? neighborlist(atoms_df(mnb.ff.system).r, unitcell=periodic_box, nonbonded_cutoff)
        : neighborlist(atoms_df(mnb.ff.system).r, nonbonded_cutoff)
    )

    vdw_interactions = Vector{BufferedVdWInteraction{T, 14, 7}}()
    es_interactions = Vector{BufferedESInteraction{T, 14, 7}}()

    vicinal_cache   = mnb.cache[:vicinal_cache]
    bond_cache      = mnb.cache[:bond_cache]
    geminal_cache   = mnb.cache[:geminal_cache]
    check_vicinal(a1, a2)   = (a1 => a2) ∈ vicinal_cache
    check_bond(a1, a2)      = (a1 => a2) ∈ bond_cache
    check_geminal(a1, a2)   = (a1 => a2) ∈ geminal_cache


    hint = Int(round(1.2 * natoms(mnb.ff.system)))
    sizehint!(vdw_interactions, hint)
    sizehint!(es_interactions, hint)

    atom_cache   = atoms(mnb.ff.system)
    idx_cache    = atoms_df(mnb.ff.system).idx
    charge_cache = atoms_df(mnb.ff.system).charge
    type_cache   = atoms_df(mnb.ff.system).atom_type

    mnb.es_constants = ESConstantTerms{T}(T(mnb.ff.options[:electrostatic_cuton]),
                T(mnb.ff.options[:electrostatic_cutoff]))

    for buf_candidate in neighbors

        buf_1 = buf_candidate[1]
        buf_2 = buf_candidate[2]

        atom_1_idx = idx_cache[buf_1]
        atom_2_idx = idx_cache[buf_2]

        #ignore direct bonds and atom pairs that are seperated by only 1 atom
        (check_bond(atom_1_idx, atom_2_idx) || check_geminal(atom_1_idx, atom_2_idx)) &&
            continue


        atom_1 = atom_cache[buf_1]
        atom_2 = atom_cache[buf_2]

        atom_1_type = type_cache[buf_1]
        atom_2_type = type_cache[buf_2]

        atom_1_charge = charge_cache[buf_1]
        atom_2_charge = charge_cache[buf_2]

        is_vicinal_pair = check_vicinal(atom_1_idx, atom_2_idx)

        rij, rij_7, eij = get_VdW_params(atom_1_type, atom_2_type, vdw_df, mnb)


        #if es_enabled
        push!(
            es_interactions,
            BufferedESInteraction{T,14,7}(
                T(buf_candidate[3]), #the distance
                es_cut_off,
                es_cut_on,
                atom_1,
                atom_2,
                atom_1_charge,
                atom_2_charge,
                is_vicinal_pair,
                ddd,
                mnb.es_switch,
                mnb.es_constants
            )
        )
        ismissing(rij) && continue
        push!(
            vdw_interactions,
            BufferedVdWInteraction{T, 14, 7}(
                T(buf_candidate[3]),
                atom_1, 
                atom_2,
                rij,
                rij_7,
                eij
            )
        )


    end

    mnb.vdw_interactions = vdw_interactions
    mnb.es_interactions = es_interactions

end 

function compute_forces(es_interaction::BufferedESInteraction{T,14,7}) where {T<:Real}

    d = es_interaction.distance
    at1 = es_interaction.a1
    at2 = es_interaction.a2
    qi = es_interaction.qi
    qj = es_interaction.qj
    is_1_4 = es_interaction.is_1_4
    direction = normalize(at1.r - at2.r)
    ddd = es_interaction.ddd
    ct = es_interaction.es_constants

    
    cut_on = es_interaction.cut_on
    cut_off = es_interaction.cut_off

    #515
    es_factor = 0
    if es_interaction.es_switch 
        d >= cut_off && return
        if !(ddd)
            if d > cut_on
                dd = (d+0.05)^2
                es_factor = ct.ac * inv(dd) + ct.bc + dd * (3* ct.c3 + 2/ct.denom * dd)
                es_factor *= ct.CONSTANT_ * qi * qj
            else
                es_factor = ct.CONSTANT_ * qi * qj *  (ddd + 1) / (1 * (d + 0.05)^(ddd + 2))
            end

        else

            if d > cut_on
                es_factor = 2 * ct.ac * inv(d + 0.05)^3 + 2 * ct.bc * inv(d + 0.05)
                es_factor += 2 *  (3 * ct.c3) * (d + 0.05) + 2 * (2/ct.denom) * (d + 0.05)^3
                es_factor *= ct.CONSTANT_ * qi * qj
            else
                es_factor = 2 * inv(d + 0.05)^3
                es_factor *= ct.CONSTANT_ * qi * qj
            end

        end

    else
        es_factor = ct.CONSTANT_ * qi * qj * (1+ddd) / (d+0.05)^(ddd + 2)
    end

    is_1_4 && (es_factor *= 0.75;)
    force = direction * es_factor * force_prefactor

    at1.F += force
    at2.F -= force
end

function compute_energy(es_interaction::BufferedESInteraction{T,14,7}) where {T<:Real}
    d = es_interaction.distance
    qi = es_interaction.qi
    qj = es_interaction.qj
    ddd = es_interaction.ddd
    is_1_4 = es_interaction.is_1_4
    ct = es_interaction.es_constants

    #here cuton and cutoff are the ES specific cutoff/cuton
    cut_on = es_interaction.cut_on
    cut_off = es_interaction.cut_off


    es = 0
    if es_interaction.es_switch 
        d >= cut_off && return 0
        dd = (d + 0.05)^2

        if !(ddd)
            if d > cut_on
                es = ct.ac - dd * (ct.bc + dd * (ct.c3 + ct.d5 * dd))
                es /= (d+ 0.05)
                es += ct.const_

                #es_factor = ct.ac * inv(r + 0.05)^2 + ct.bc + dd * (3*ct.c3)
            else
                es = inv(d + 0.05) + ct.add
            end
        else
            if d > cut_on
                es =  (ct.ac * inv(dd) + 2 * ct.bc * log(inv(d+0.05)) 
                    - dd * (3*ct.c3 + dd/ct.denom) + ct.constr_)
            else
                es = inv(dd) + ct.addr
            end
        end

        es *= ct.CONSTANT_ * qi * qj
    else
        es = ct.CONSTANT_ * qi * qj / (d + 0.05)^(1 + ddd)
    end

    is_1_4 && (es *= 0.75;)
    return es

end

function compute_energy(vdw_interaction::BufferedVdWInteraction{T,14,7}) where {T<:Real}
    vdw_interaction.distance == 0 && return     #is this needed with cellListMap.jl?
    d     = vdw_interaction.distance
    eij   = vdw_interaction.eij
    rij   = vdw_interaction.rij
    rij_7 = vdw_interaction.rij_7

    #MMFF94 Paper2: Eq1
    left = eij * ((1.07*rij) / (d + 0.07 * rij))^7
    right = ((1.12 * rij_7) / (d^7 + 0.12 * rij_7)) - 2

    e = (left * right)u"cal" |> u"J"

    ustrip(e)
end 


function compute_forces(vdw_interaction::BufferedVdWInteraction{T,14,7}) where {T<:Real}
    #i assume this is the actual distance, i.e. the square root of the substraction    d = vdw_interaction.distance
    at1     = vdw_interaction.a1
    at2     = vdw_interaction.a2
    eij     = vdw_interaction.eij
    rij     = vdw_interaction.rij
    d       = vdw_interaction.distance
    direction = normalize(at1.r - at2.r)

    dbuf = 0.07
    gbuf = 0.12
    

    q = d / rij
    p = (1 + dbuf) / (q + dbuf)
    h = (1 + gbuf) / (q^7 + gbuf)
    gbuf = 0.12
    dbuf = 0.07
    fac   = (eij * p^7 * (-7.))
    left  = (q^6 * h) / (q^7 + gbuf)
    right = ((h-2.)) / (q+dbuf)

    vdw_factor = ((fac * (left + right)) / rij) 
    #this is according to C++ BALL's implementation
    force = direction * vdw_factor * force_prefactor * joule_per_cal


    at1.F -= force
    at2.F += force

end

function compute_forces(mnb::MNonBondedComponent{T}) where {T<:Real}
    constrained_ids = getproperty.(atoms(mnb.ff.system)[mnb.ff.constrained_atoms], :idx)

    filter_pairs = (
        isempty(constrained_ids) 
            ? identity 
            : s -> filter(p -> (p.a1.idx ∉ constrained_ids) || (p.a2.idx ∉ constrained_ids), s)
    )

    #compute_forces(mnb.vdw_interactions[1])
    get(mnb.ff.options ,:MMFF_VDW_ENABLED, true) && map(compute_forces, filter_pairs(mnb.vdw_interactions))
    get(mnb.ff.options ,:MMFF_ES_ENABLED , true) && map(compute_forces, filter_pairs(mnb.es_interactions))

    nothing
end

function compute_energy(mnb::MNonBondedComponent{T}) where {T<:Real}
    # iterate over all interactions in the system
    compute_energy(mnb.es_interactions[1])
    vdw_energy   = get(mnb.ff.options ,:MMFF_VDW_ENABLED, true) ? mapreduce(compute_energy, +, mnb.vdw_interactions; init=zero(T)) : 0
    es_energy    = get(mnb.ff.options ,:MMFF_ES_ENABLED , true) ? mapreduce(compute_energy, +, mnb.es_interactions;  init=zero(T)) : 0

    mnb.energy["Van der Waals"]  = vdw_energy
    mnb.energy["Electrostatic"]  = es_energy
    
    vdw_energy + es_energy
end