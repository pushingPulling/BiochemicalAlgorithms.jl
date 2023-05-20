export LennardJonesInteraction, ElectrostaticInteraction, NonBondedComponent


# e_scaling_factor contains the unit conversions und the constants
# appearing in Coulomb's law:
#
#               1        q1 * e0 * q2 * e0
#     F = ------------- ------------------
#         4 PI epsilon0       r * r
#
# Conversion factors are 1e-10 for Angstrom -> m
# and Constants::e0 for the proton charge
const ES_Prefactor = ustrip(Constants.N_A * Constants.e₀^2 / (4.0 * π * Constants.ε₀) |> u"angstrom*kJ/mol")

const ES_Prefactor_force = ustrip(Constants.e₀^2 / (4.0 * π * Constants.ε₀) * 1u"angstrom^-2" |> u"N")

@auto_hash_equals struct CubicSwitchingFunction{T<:Real}
    cutoff::T
    cuton::T

    sq_cutoff::T
    sq_cuton::T

    inverse_distance_off_on_3::T

    function CubicSwitchingFunction{T}(cutoff::T, cuton::T) where {T<:Real}
        sq_cutoff = cutoff^2
        sq_cuton  = cuton^2
    
        inverse_distance_off_on_3 = (sq_cutoff - sq_cuton)^3
    
        if (inverse_distance_off_on_3 <= T(0.0))
            @warn "NonBondedComponent(): cuton value should be smaller than cutoff -- " *
                    "switching function disabled."
    
            inverse_distance_off_on_3 = T(0.0)
        else
            inverse_distance_off_on_3 = T(1.0) / inverse_distance_off_on_3
        end

        new(cutoff, cuton, sq_cutoff, sq_cuton, inverse_distance_off_on_3)
    end
end

# the switching function is defined as follows:
#         (r_{off}^2 - R^2)^2 (r_{off}^2 + 2 R^2 - 3r_{on}^2)
# sw(R) = ---------------------------------------------------
#                    (r_{off}^2 - r_{on}^2)^3
#
# [Brooks et al., J. Comput. Chem., 4:191 (1983)]
#
# the derivative has the following form:
#
#                      (r_{off}^2 - R^2)(r_{on}^2 - R^2)
#   d/dR sw(R) = 12 R -----------------------------------
#                           (r_{off}^2 - r_{on}^2)^3
#
function (f::CubicSwitchingFunction{T})(sq_distance::T) where {T<:Real}
    below_off = ((sq_distance < f.sq_cutoff) ? one(T) : zero(T));
    below_on  = ((sq_distance < f.sq_cuton)  ? one(T) : zero(T));

    below_off * (below_on + (one(T) - below_on) * (f.sq_cutoff - sq_distance)^2
                * (f.sq_cutoff + T(2.0) * sq_distance - T(3.0) * f.sq_cuton)
                * f.inverse_distance_off_on_3)
end

function switching_derivative(f::CubicSwitchingFunction{T}, sq_distance::T) where {T<:Real}
    primal = (f.sq_cutoff - sq_distance)^2 *
                (f.sq_cutoff + T(2.0) * sq_distance - T(3.0) * f.sq_cuton) *
                f.inverse_distance_off_on_3

    difference_to_off = f.sq_cutoff - sq_distance
    difference_to_on  = f.sq_cuton  - sq_distance

    derivative = T(12.0) * difference_to_off * difference_to_on *
                    f.inverse_distance_off_on_3

    primal, derivative
end

@auto_hash_equals struct LennardJonesInteraction{T<:Real, N, M}
    A::T
    B::T
    distance::T
    scaling_factor::T
    a1::Atom{T}
    a2::Atom{T}
    switching_function
end

@auto_hash_equals struct ElecrostaticInteraction{T<:Real}
    q1q2::T
    distance::T
    scaling_factor::T
    distance_dependent_dielectric::Bool
    a1::Atom{T}
    a2::Atom{T}
    switching_function
end

@inline function _try_assign_vdw!(
        atom_1::Atom{T},
        type_atom_1,
        atom_2::Atom{T},
        type_atom_2,
        distance::T,
        scaling_factor::T,
        lj_combinations,
        lj_interactions,
        switching_function,
        ff::ForceField{T}) where {T<:Real}

    params = get(lj_combinations, (I=type_atom_1, J=type_atom_2,), missing)

    if !ismissing(params)
        push!(
            lj_interactions,
                LennardJonesInteraction{T, 12, 6}(
                    params.A_ij,
                    params.B_ij,
                    T(distance),
                    scaling_factor,
                    atom_1,
                    atom_2,
                    switching_function
                )
        )
    else
        @warn "NonBondedComponent(): cannot find vdW parameters for "                   *
                "atom types $(type_atom_1)-$(type_atom_2) (atoms are: "                 *
                "$(get_full_name(atom_1, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID))/" *
                "$(get_full_name(atom_2, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID)))"

        !(atom_1 in ff.unassigned_atoms) && push!(ff.unassigned_atoms, atom_1)
        !(atom_2 in ff.unassigned_atoms) && push!(ff.unassigned_atoms, atom_2)

        if length(ff.unassigned_atoms) > ff.options[:max_number_of_unassigned_atoms]
            throw(TooManyErrors())
        end
    end
end

@auto_hash_equals mutable struct NonBondedComponent{T<:Real} <: AbstractForceFieldComponent{T}
    name::String
    ff::ForceField{T}
    cache::Dict{Symbol, Any}
    energy::Dict{String, T}
    lj_interactions::AbstractVector{LennardJonesInteraction{T, 12, 6}}
    hydrogen_bonds::AbstractVector{LennardJonesInteraction{T, 12, 10}}
    electrostatic_interactions::AbstractVector{ElecrostaticInteraction{T}}

    function NonBondedComponent{T}(ff::ForceField{T}) where {T<:Real}
        this = new("NonBonded", ff, Dict{Symbol, Any}(), Dict{String, T}())

        setup!(this)
        update!(this)

        this
    end
end

function _setup_vdw!(nbc::NonBondedComponent{T}) where {T<:Real}
    # extract the parameter section for quadratic angle bends
    lj_section = extract_section(nbc.ff.parameters, "LennardJones")
    lj_df = lj_section.data

    unit_R = get(lj_section.properties, "unit_R",       "angstrom")
    unit_ϵ = get(lj_section.properties, "unit_epsilon", "kcal/mol")

    # ball used to write Angstrom with a capital letter; this clashes with the convention in Unitful.jl
    if strip(unit_R) == "Angstrom"
        unit_R = "angstrom"
    end

    R_factor = ustrip((1uparse(unit_R)) |> u"angstrom")
    ϵ_factor = ustrip((1uparse(unit_ϵ)) |> u"kJ/mol")

    lj_minimal_df = lj_df[:, [:I, :R, :epsilon]]
    lj_minimal_df.R *= R_factor
    lj_minimal_df.epsilon *= ϵ_factor

    lj_combinations = crossjoin(lj_minimal_df, lj_minimal_df; makeunique=true)
    rename!(lj_combinations, 
        [:R => :R_i, :epsilon => :epsilon_i, :I_1 => :J, :R_1 => :R_j, :epsilon_1 => :epsilon_j])

    # convert the data from RE - format to A/B - format
    radius_averaging = get(lj_section.properties, :radius_averaging, "arithmetic")

    if radius_averaging ∉ ["arithmetic", "geometric"]
        @warn "NonBondedComponent(): unknown radius averaging method $(radius_averaging)."
        @warn "Using arithmetic averaging instead."

        radius_averaging = "arithmetic"
    end

    r_6 = ((radius_averaging == "arithmetic")
            ? (lj_combinations.R_i .+ lj_combinations.R_j)
            : T(2.0) * sqrt.(lj_combinations.R_i .* lj_combinations.Rj)
        ).^6

    ϵ = sqrt.(lj_combinations.epsilon_i .* lj_combinations.epsilon_j)

    lj_combinations.A_ij = ϵ .* r_6.^2
    lj_combinations.B_ij = T(2.0) * ϵ .* r_6
    
    lj_combinations_cache = Dict(
        (I=r.I, J=r.J,) => (A_ij=r.A_ij, B_ij=r.B_ij) 
            for r in eachrow(lj_combinations)
    )

    vdw_cutoff       = nbc.ff.options[:vdw_cutoff]
    vdw_cuton        = nbc.ff.options[:vdw_cuton]

    # build the switching function as a closure
    vdw_switching_function = CubicSwitchingFunction{T}(vdw_cutoff, vdw_cuton)
  
    scaling_vdw_1_4 = nbc.ff.options[:scaling_vdw_1_4]
    if (scaling_vdw_1_4 == T(0.0))
        @warn "NonBondedComponent(): illegal - 1-4 vdW scaling factor: must be non-zero!"
        @warn "Resetting to 1.0."

        scaling_vdw_1_4 = T(1.0)
    else
        scaling_vdw_1_4 = T(1.0) / scaling_vdw_1_4
    end

    # remember those parts that stay constant when only the system is updated
    nbc.cache[:lj_combinations] = lj_combinations_cache

    nbc.cache[:vdw_cutoff] = vdw_cutoff
    nbc.cache[:vdw_cuton]  = vdw_cuton

    nbc.cache[:vdw_switching_function] = vdw_switching_function

    nbc.cache[:scaling_vdw_1_4] = scaling_vdw_1_4
end

function _setup_hydrogenbonds!(nbc::NonBondedComponent{T}) where {T<:Real}
    # we will also need the hydrogen bond section to determine if parameters are missing for a pair of atoms
    hydrogen_bonds_section = extract_section(nbc.ff.parameters, "HydrogenBonds")
    hydrogen_bonds_df = hydrogen_bonds_section.data

    unit_hbond_A = get(hydrogen_bonds_section.properties, "unit_A", "kcal/mol*angstrom^12")
    unit_hbond_B = get(hydrogen_bonds_section.properties, "unit_B", "kcal/mol*angstrom^10")

    # here, BALL used A for angstrom; this clashes with the convention in Unitful.jl
    unit_hbond_A = replace(unit_hbond_A, "A" => "angstrom")
    unit_hbond_B = replace(unit_hbond_B, "A" => "angstrom")

    hbond_A_factor = ustrip((1uparse(unit_hbond_A)) |> u"kJ/mol*angstrom^12")
    hbond_B_factor = ustrip((1uparse(unit_hbond_B)) |> u"kJ/mol*angstrom^10")

    # hbond terms are already stored in A/B - format, so conversion is simple
    hydrogen_bonds_df.A .*= T(hbond_A_factor)
    hydrogen_bonds_df.B .*= T(hbond_B_factor)

    hydrogen_bond_combinations_cache = Dict(
        (I=r.I, J=r.J,) => (A=r.A, B_ij=r.B) 
            for r in eachrow(hydrogen_bonds_df)
    )

    # remember those parts that stay constant when only the system is updated
    nbc.cache[:hydrogen_bond_combinations] = hydrogen_bond_combinations_cache
end

function _setup_electrostatic_interactions!(nbc::NonBondedComponent{T}) where {T<:Real}
    es_cutoff        = nbc.ff.options[:electrostatic_cutoff]
    es_cuton         = nbc.ff.options[:electrostatic_cuton]  

    es_switching_function  = CubicSwitchingFunction{T}(es_cutoff, es_cuton)
   
    scaling_es_1_4 = nbc.ff.options[:scaling_electrostatic_1_4]
    if (scaling_es_1_4 == T(0.0))
        @warn "NonBondedComponent(): illegal - 1-4 electrostatic scaling factor: must be non-zero!"
        @warn "Resetting to 1.0."

        scaling_es_1_4 = T(1.0)
    else
        scaling_es_1_4 = T(1.0) / scaling_es_1_4
    end

    # remember those parts that stay constant when only the system is updated
    nbc.cache[:es_cutoff] = es_cutoff
    nbc.cache[:es_cuton]  = es_cuton

    nbc.cache[:es_switching_function] = es_switching_function

    nbc.cache[:scaling_es_1_4] = scaling_es_1_4
end

function setup!(nbc::NonBondedComponent{T}) where {T<:Real}
    _setup_vdw!(nbc)
    _setup_hydrogenbonds!(nbc)
    _setup_electrostatic_interactions!(nbc)

    # the cutoffs for the nonbonded pair list and the switching function
    nonbonded_cutoff = nbc.ff.options[:nonbonded_cutoff]

    # do we need to construct a periodic box?
    periodic_box = [
        nbc.ff.options[:periodic_box_width], 
        nbc.ff.options[:periodic_box_height],
        nbc.ff.options[:periodic_box_depth]
    ]

    # cache the is_bound_to-, geminal-, and vicinal-relationships
    mg = MetaGraph(non_hydrogen_bonds_df(nbc.ff.system), :a1, :a2)

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
    nbc.cache[:nonbonded_cutoff] = nonbonded_cutoff
    nbc.cache[:periodic_box]     = periodic_box

    nbc.cache[:bond_cache]    = bond_cache
    nbc.cache[:geminal_cache] = geminal_cache
    nbc.cache[:vicinal_cache] = vicinal_cache
end

function update!(nbc::NonBondedComponent{T}) where {T<:Real}
    ff = nbc.ff

    periodic_box = nbc.cache[:periodic_box]
    nonbonded_cutoff = nbc.cache[:nonbonded_cutoff]

    scaling_vdw_1_4 = nbc.cache[:scaling_vdw_1_4]
    scaling_es_1_4  = nbc.cache[:scaling_es_1_4]

    vdw_switching_function = nbc.cache[:vdw_switching_function]
    es_switching_function  = nbc.cache[:es_switching_function]

    lj_combinations = nbc.cache[:lj_combinations]

    hydrogen_bond_combinations = nbc.cache[:hydrogen_bond_combinations]

    bond_cache    = nbc.cache[:bond_cache]
    geminal_cache = nbc.cache[:geminal_cache]
    vicinal_cache = nbc.cache[:vicinal_cache]

    check_bond(a1, a2)    = (a1 => a2) ∈ bond_cache
    check_geminal(a1, a2) = (a1 => a2) ∈ geminal_cache
    check_vicinal(a1, a2) = (a1 => a2) ∈ vicinal_cache

    neighbors = ((ff.options[:periodic_boundary_conditions]) 
        ? neighborlist(atoms_df(ff.system).r, unitcell=periodic_box, nonbonded_cutoff)
        : neighborlist(atoms_df(ff.system).r, nonbonded_cutoff)
    )

    distance_dependent_dielectric = ff.options[:distance_dependent_dielectric]

    lj_interactions = Vector{LennardJonesInteraction{T, 12, 6}}()
    hydrogen_bonds  = Vector{LennardJonesInteraction{T, 12, 10}}()
    electrostatic_interactions = Vector{ElecrostaticInteraction{T}}()

    hint = Int(round(1.2 * natoms(ff.system)))
    sizehint!(lj_interactions, hint)
    sizehint!(electrostatic_interactions, hint)

    atom_cache   = atoms(ff.system)
    idx_cache    = atoms_df(ff.system).idx
    charge_cache = atoms_df(ff.system).charge
    type_cache   = atoms_df(ff.system).atom_type

    for lj_candidate in neighbors
        lj_1 = lj_candidate[1]
        lj_2 = lj_candidate[2]

        atom_1 = atom_cache[lj_1]
        atom_2 = atom_cache[lj_2]

        atom_1_idx = idx_cache[lj_1]
        atom_2_idx = idx_cache[lj_2]

        atom_1_type = type_cache[lj_1]
        atom_2_type = type_cache[lj_2]

        # exclude 1-2 and 1-3 interactions
        if check_bond(atom_1_idx, atom_2_idx) || check_geminal(atom_1_idx, atom_2_idx)
            continue
        end

        vicinal_pair = check_vicinal(atom_1_idx, atom_2_idx)

        q1q2 = charge_cache[lj_1]*charge_cache[lj_2]

        if q1q2 ≠ zero(T)
            push!(
                electrostatic_interactions,
                ElecrostaticInteraction{T}(
                    q1q2,
                    T(lj_candidate[3]),
                    vicinal_pair ? scaling_es_1_4 : T(1.0),
                    distance_dependent_dielectric,
                    atom_1,
                    atom_2,
                    es_switching_function
                )
            )
        end

        # first, figure out if the atoms are part of a torsion
        if !vicinal_pair              
            # no. now, figure out if these atoms are part of a hydrogen bond

            # parameters for hydrogen bonds are used, if they exist
            # and the two atoms are not vicinal (1-4).

            h_params = coalesce(
                get(hydrogen_bond_combinations, (I=atom_1_type, J=atom_2_type,), missing),
                get(hydrogen_bond_combinations, (I=atom_2_type, J=atom_1_type,), missing)
            )

            if !ismissing(h_params)
                params = only(params)
                push!(
                    hydrogen_bonds,
                    LennardJonesInteraction{T, 12, 10}(
                        params.A,
                        params.B,
                        T(lj_candidate[3]),
                        T(1.0),
                        atom_1,
                        atom_2,
                        vdw_switching_function
                    )
                )
            else
                _try_assign_vdw!(
                    atom_1,
                    atom_1_type,
                    atom_2,
                    atom_2_type,
                    T(lj_candidate[3]),
                    T(1.0),
                    lj_combinations,
                    lj_interactions,
                    vdw_switching_function,
                    nbc.ff
                )
            end
        else
            # this is a torsion
            _try_assign_vdw!(
                atom_1,
                atom_1_type,
                atom_2,
                atom_2_type,
                T(lj_candidate[3]),
                scaling_vdw_1_4,
                lj_combinations,
                lj_interactions,
                vdw_switching_function,
                nbc.ff
            )
        end
    end

    nbc.lj_interactions            = lj_interactions
    nbc.hydrogen_bonds             = hydrogen_bonds
    nbc.electrostatic_interactions = electrostatic_interactions

    nothing
end

@inline function compute_energy(lji::LennardJonesInteraction{T, 12, 6}) where {T<:Real}
    inv_dist_6 = lji.distance^-6

    inv_dist_6 * (inv_dist_6 * lji.A - lji.B) * lji.scaling_factor * lji.switching_function(lji.distance^2)
end

@inline function compute_energy(hb::LennardJonesInteraction{T, 12, 10}) where {T<:Real}
    hb.distance^-12 * hb.A - hb.distance^-10 * hb.B * hb.scaling_factor * hb.switching_function(hb.distance^2)
end

@inline function compute_energy(esi::ElecrostaticInteraction{T}) where {T<:Real}
    energy = esi.distance_dependent_dielectric ? 0.25 * esi.q1q2 / (esi.distance^2) : esi.q1q2 / esi.distance

    energy * esi.scaling_factor * esi.switching_function(esi.distance^2) * T(ES_Prefactor)
end

function compute_energy(nbc::NonBondedComponent{T}) where {T<:Real}
    # iterate over all interactions in the system
    vdw_energy   = mapreduce(compute_energy, +, nbc.lj_interactions;            init=zero(T))
    hbond_energy = mapreduce(compute_energy, +, nbc.hydrogen_bonds;             init=zero(T))
    es_energy    = mapreduce(compute_energy, +, nbc.electrostatic_interactions; init=zero(T))

    nbc.energy["Van der Waals"]  = vdw_energy
    nbc.energy["Hydrogen Bonds"] = hbond_energy
    nbc.energy["Electrostatic"]  = es_energy
    
    vdw_energy + hbond_energy + es_energy
end

function compute_forces(lji::LennardJonesInteraction{T, 12, 6}) where {T<:Real}
    direction = lji.a1.r - lji.a2.r

    sq_distance = squared_norm(direction)

    if (sq_distance > zero(T) && sq_distance <= lji.switching_function.sq_cutoff)
        factor = T(force_prefactor) / sq_distance
        inv_distance_6 = sq_distance^-3

        factor *= inv_distance_6 * lji.scaling_factor * (12 * lji.A * inv_distance_6 - 6 * lji.B)

        # do we have to use the switching function (cuton <= distance <= cutoff)?
        if (sq_distance > lji.switching_function.sq_cuton)
            switch_value, switch_derivative = switching_derivative(lji.switching_function, sq_distance)

            # First, multiply the current force with the switching function
            factor *= switch_value

            # Second, we add the product of the energy and the derivative
            # of the switching function (the total force is the derivative of
            # a product of functions)
            energy = -T(force_prefactor) * lji.scaling_factor *
                inv_distance_6 * (inv_distance_6 * lji.A - lji.B)
  
            factor += switch_derivative * energy
        end
    
        force = factor * direction

        #@info "vdW $(get_full_name(lji.a1))<->$(get_full_name(lji.a2)) $(force)"

        lji.a1.F += force
        lji.a2.F -= force
    end
end

function compute_forces(hb::LennardJonesInteraction{T, 12, 10}) where {T<:Real}
    direction = hb.a1.r - hb.a2.r

    sq_distance = squared_norm(direction)

    if (sq_distance > zero(T) && sq_distance <= hb.switching_function.sq_cutoff)
        inv_distance_2 =  T(1.0) / sq_distance
        inv_distance_10 = sq_distance^-5
        inv_distance_12 = sq_distance^-6
      
        factor = T(force_prefactor) * inv_distance_2
        
        factor *= inv_distance_12 * (12 * hb.A * inv_distance_2 - 10 * hb.B); 

        # do we have to use the switching function (cuton <= distance <= cutoff)?
        if (sq_distance > hb.switching_function.sq_cuton)
            switch_value, switch_derivative = switching_derivative(hb.switching_function, sq_distance)

            # First, multiply the current force with the switching function
            factor *= switch_value

            # Second, we add the product of the energy and the derivative
            # of the switching function (the total force is the derivative of
            # a product of functions)
            energy = -T(force_prefactor) * hb.scaling_factor *
                inv_distance_10 * (hb.A * inv_distance_2 - hb.B)
  
            factor += switch_derivative * energy
        end
    
        force = factor * direction

        @info "HB $(get_full_name(hb.a1))<->$(get_full_name(hb.a2)) $(force)"


        hb.a1.F += force
        hb.a2.F -= force
    end
end

function compute_forces(esi::ElecrostaticInteraction{T}) where {T<:Real}
    direction = esi.a1.r - esi.a2.r

    sq_distance = squared_norm(direction)

    if (sq_distance > zero(T) && sq_distance <= esi.switching_function.sq_cutoff)
        inv_distance_2 = T(1.0) / sq_distance
        inv_distance = sqrt(inv_distance_2)
      
        factor = esi.q1q2 * inv_distance_2 * esi.scaling_factor * ES_Prefactor_force

        # distinguish between constant and distance dependent dielectric 
        if esi.distance_dependent_dielectric
            # distance dependent dielectric:  epsilon = 4 * r_ij
            # 4 reduces to 2 (due to derivation of the energy)
            factor *= 0.5 * inv_distance_2
        else 
            # distance independent dielectric constant
            factor *= sqrt(inv_distance_2)
        end

        # we have to use the switching function (cuton <= distance <= cutoff)
        if (sq_distance > esi.switching_function.sq_cuton)
            switch_value, switch_derivative = switching_derivative(esi.switching_function, sq_distance)

            # First, multiply the current force with the switching function
            factor *= switch_value

            # Second, we add the product of the energy and the derivative
            # of the switching function (the total force is the derivative of
            # a product of functions)
            # we save the multiplication by distance to avoid the normalization
            # of the direction vector
            # In fact, we calculate the negative energy, since we de not
            # calculate the force, but the derivative of the energy above.

            # calculate the electrostatic energy
            dist_depend_factor = ((esi.distance_dependent_dielectric) 
                ? T(0.25) * inv_distance : one(T))
            
            energy = -T(ES_Prefactor_force) * esi.scaling_factor * dist_depend_factor *
                inv_distance * esi.q1q2

            factor += switch_derivative * energy
        end
    
        force = factor * direction

        # @info "ES force $(get_full_name(esi.a1))<->$(get_full_name(esi.a2)) $(force)"
        
        esi.a1.F += force
        esi.a2.F -= force
    end
end

function compute_forces(nbc::NonBondedComponent{T}) where {T<:Real}
    constrained_ids = getproperty.(atoms(nbc.ff.system)[nbc.ff.constrained_atoms], :idx)

    filter_pairs = (
        isempty(constrained_ids) 
            ? identity 
            : s -> filter(p -> (p.a1.idx ∉ constrained_ids) || (p.a2.idx ∉ constrained_ids), s)
    )

    map(compute_forces, filter_pairs(nbc.lj_interactions))
    map(compute_forces, filter_pairs(nbc.hydrogen_bonds))
    map(compute_forces, filter_pairs(nbc.electrostatic_interactions))

    nothing
end
