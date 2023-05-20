export QuadraticBondStretch, QuadraticStretchComponent

@auto_hash_equals struct QuadraticBondStretch{T<:Real}
    r0::T
    k::T
    a1::Atom{T}
    a2::Atom{T}
end

@auto_hash_equals mutable struct QuadraticStretchComponent{T<:Real} <: AbstractForceFieldComponent{T}
    name::String
    ff::ForceField{T}
    cache::Dict{Symbol, Any}
    energy::Dict{String, T}
    stretches::AbstractVector{QuadraticBondStretch{T}}

    function QuadraticStretchComponent{T}(ff::ForceField{T}) where {T<:Real}
        this = new("QuadraticBondStretch", ff, Dict{Symbol, Any}(), Dict{String, T}())

        setup!(this)
        update!(this)

        this
    end
end

function setup!(qsc::QuadraticStretchComponent{T}) where {T}
    ff = qsc.ff
    
    # extract the parameter section for quadratic bond stretches
    stretch_section = extract_section(qsc.ff.parameters, "QuadraticBondStretch")
    stretch_df = stretch_section.data

    unit_k  = get(stretch_section.properties, "unit_k", "kcal/mol")
    unit_r0 = get(stretch_section.properties, "unit_r0", "angstrom")

    # ball used to write Angstrom with a capital letter; this clashes with the convention in Unitful.jl
    if strip(unit_r0) == "Angstrom"
        unit_r0 = "angstrom"
    end

    k_factor  = ustrip((1uparse(unit_k))  |> u"kJ/mol")
    r0_factor = ustrip((1uparse(unit_r0)) |> u"angstrom")

    # group the stretch parameters by type_i, type_j combinations
    stretch_combinations = groupby(stretch_df, ["I", "J"])

    stretches = Vector{QuadraticBondStretch}(undef, nbonds(ff.system))

    # iterate over all bonds in the system
    for (i, bond) in enumerate(bonds(ff.system))
        a1 = atom_by_idx(parent_system(ff.system), bond.a1)
        a2 = atom_by_idx(parent_system(ff.system), bond.a2)

        type_a1 = a1.atom_type
        type_a2 = a2.atom_type

        if has_flag(bond, :TYPE__HYDROGEN)
            # skip hydrogen bonds
            stretches[i] = QuadraticBondStretch(one(T), zero(T), a1, a2)
        end

        qbs = coalesce(
            get(stretch_combinations, (I=type_a1, J=type_a2,), missing),
            get(stretch_combinations, (I=type_a2, J=type_a1,), missing),
            get(stretch_combinations, (I=type_a1, J="*",    ), missing),
            get(stretch_combinations, (I="*",     J=type_a2,), missing),
            get(stretch_combinations, (I="*",     J="*",    ), missing)
        )

        stretches[i] = if ismissing(qbs)

            @warn "QuadraticStretchComponent(): cannot find stretch parameters for " *
                "atom types $(type_a1)-$(type_a2) (atoms are: "                      *
                "$(get_full_name(a1, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID))/"  *
                "$(get_full_name(a2, FullNameType.ADD_VARIANT_EXTENSIONS_AND_ID)))"

            push!(ff.unassigned_atoms, a1)
            push!(ff.unassigned_atoms, a2)

            if length(ff.unassigned_atoms) > ff.options[:max_number_of_unassigned_atoms]
                throw(TooManyErrors())
            end

            # we don't want to get any force or energy component from this stretch
            QuadraticBondStretch(one(T), zero(T), a1, a2)
        else
            QuadraticBondStretch(
                T(r0_factor*only(qbs.r0)), 
                T(k_factor*only(qbs.k)),
                a1,
                a2
            )
        end
    end

    qsc.stretches = stretches
end

function update!(qsc::QuadraticStretchComponent{T}) where {T<:Real}  
    nothing
end

@inline function compute_energy(qbs::QuadraticBondStretch{T}) where {T<:Real}
    d = distance(qbs.a1.r, qbs.a2.r)

    qbs.k * (d - qbs.r0)^2
end

function compute_energy(qsc::QuadraticStretchComponent{T}) where {T<:Real}
    # iterate over all bonds in the system

    total_energy = mapreduce(compute_energy, +, qsc.stretches; init=zero(T))

    qsc.energy["Bond Stretches"] = total_energy

    total_energy
end

function compute_forces(qbs::QuadraticBondStretch{T}) where {T<:Real}
    direction = qbs.a1.r - qbs.a2.r
    distance = norm(direction)

    if distance == zero(T)
        return
    end

    direction *= force_prefactor * 2 * qbs.k * (distance - qbs.r0) / distance

    #@info "$(get_full_name(qbs.a1))<->$(get_full_name(qbs.a2)) $(direction)"
    
    qbs.a1.F -= direction
    qbs.a2.F += direction
end

function compute_forces(qsc::QuadraticStretchComponent{T}) where {T<:Real}
    map(compute_forces, qsc.stretches)

    nothing
end