@auto_hash_equals mutable struct MBend{T} <: AbstractForceFieldComponent{T}
    theta0   ::T
    thata_d  ::T
    theta    ::T
    ka       ::T
    at1      ::Atom{T}
    at2      ::Atom{T}
    at3      ::Atom{T}
    is_linear::Bool
    ATIJK    ::Int
    n1       ::Vector3{T}
    n2       ::Vector3{T}
end

@auto_hash_equals mutable struct MStretch{T} <: AbstractForceFieldComponent{T}
    at1      ::Atom{T}
    at2      ::Atom{T}
    kb       ::T
    r0       ::T
    delta_r  ::T
    sbmb     ::Bool
    empirical::Bool
    #n        ::Vector3{T}
end

@auto_hash_equals mutable struct MStretchBend{T} <: AbstractForceFieldComponent{T}
    kba_ijk     ::T
    kba_kji     ::T
    sbtijk      ::T
    stretch_i_j ::MStretch
    stretch_j_k ::MStretch
    bend        ::MBend
end

@auto_hash_equals mutable struct MStretchBendComponent{T} <: AbstractForceFieldComponent{T}
    ff              ::ForceField{T}
    name            ::String
    cache           ::Dict{Symbol, Any}
    energy          ::Dict{String, T}
    # the atoms are sorted by atom_idx : lowest first
    stretches       ::Dict{NTuple{2, Atom{T}}, MStretch{T}}
    bends           ::AbstractVector{MBend}
    stretch_bends   ::AbstractVector{MStretchBend}    
end

function MStretchBendComponent{T}(ff)
    MStretchBendComponent{Float32}(f, "msbc", Dict{Symbol, Any}(),
        Dict{String, T}(), Dict{NTuple{2,Atom{T}, MStretch{T}}}(), 
        MBend{T}[], MStretchBend{T}[])
end

struct Stretch_Idx
    b::Int8
    c::Int8

    function Stretch_Idx(b::AbstractString, c::AbstractString)
        return b < c ? new(parse(Int8,b),parse(Int8,c)) : new(parse(Int8,c),parse(Int8,b))
    end

    function Stretch_Idx(b::Int, c::Int)
        return b < c ? new(b,c) : new(c,b)
    end
end

SI(b,c) = Stretch_Idx(b,c)
at_idx(type) = parse(Int8, type) > 82 ? parse(Int8, type) - 4 : parse(Int8, type)

function setup!(msb::MStretchBendComponent{T}, aro_rings, all_rings) where T<:Real
    unassigned_atoms   = Atom{T}[]
    push!(unassigned_atoms, setupStretches(msb, aro_rings)...)
    @show "stretches donee"
    setupBends(msb, all_rings)
    @show "bend donee"
    setupStretchBends(msb)
    @show "finished setup"
end

function calculate_stretch_r0(bond, at1, at2, radii,
            electro_neg, atom_type_properties, aro_rings)

    # types 83-86 dont exist in the atoms type file.
    # this func allows indexing the atom type properties df
    

    a1_n = Int(at1.element)
    a2_n = Int(at2.element)

    # default is same value as in mmff94_parameters.jl::assign_to
    at1.atom_type == "0" || at2.atom_type == "0" && return -1.0
    a1_n > 53 || a2_n > 53 || a1_n == 0 || a2_n == 0 && return -1.0

    r1 = radii[a1_n]
    r2 = radii[a2_n]
    r1 == 0.0 || r2 == 0.0 && return -1.0

    if a1_n != 0 && a1_n != 0
        bond.order in (BondOrder.Unknown, BondOrder.Quadruple) && return -1.0
        b1 = atom_type_properties[at_idx(at1.atom_type), :mltb]
        b2 = atom_type_properties[at_idx(at2.atom_type), :mltb]

        if     b1 == b2 == 0;    bond.order = BondOrder.Quadruple;
        elseif b1  + b2 == 3;    bond.order = BondOrder.Aromatic;
        else
            if any(r -> at1 in r && at2 in r, aro_rings)
                !Bool(atom_type_properties[at_idx(at1.type), :pilp]) &&
                    !Bool(atom_type_properties[at_idx(at2.type), :pilp]) ?
                        bond.order = 4 : bond.order = 5
            end
        end

        if     bond.order == 5; r1 -= 0.04 ;    r2 -= 0.04 ;
        elseif bond.order == 4; r1 -= 0.075;    r2 -= 0.075;
        elseif bond.order == 3; r1 -= 0.17 ;    r2 -= 0.17 ;
        elseif bond.order == 2; r1 -= 0.1  ;    r2 -= 0.1  ;
        else   #bond.order == 0

            h1 = 3
            h2 = 3
            if      b1 == 1 || b1 == 2; h1 = 2
            elseif  b1 == 3;            h1 = 1 end

            if      b2 == 1 || b2 == 2; h2 = 2
            elseif  b2 == 3;            h1 = 1 end

            if      h1 == 1; r1 -= 0.08
            elseif  h1 == 2; r1 -= 0.03 end

            if      h2 == 1; r2 -= 0.08
            elseif  h2 == 2; r2 -= 0.03 end

        end

    end

    #calculate shrink factor
    d = 0.008

    (a1_n == 1 || a2_n == 1) && (d = 0.0)
    (a1_n > 10 || a2_n > 10) && (d = 0.0)

    # c and n are constants defined in R.Blom and A. Haaland,
    # J. Molec. Struc, 1985, 128, 21-27.
	# calculate proportionality constant c
    (a1_n == 1 || a2_n == 1) ? c = 0.05 : c = 0.08
    n = 1.4

    delta_e = abs(electro_neg[a1_n] - electro_neg[a2_n])
    r0 = r1 + r2 - (c * n^delta_e - d)
    return r0 

end


function calculate_stretch_kb(r0, at1, at2, gdf_stretch_e)
    row = get(gdf_stretch_e, (type1=Int8(at1.element), type2=Int8(at1.element)), missing)
    if !ismissing(row) 
        kb = only(row.kb) * (only(row.r0)/r0)^6 
        return kb
    end

    at1_n = Int8(at1.element)
    at2_n = Int8(at2.element)
    at1_n < at2_n ? (p1 = at1_n; p2 = at2_n) : (p2 = at1_n; p1 = at2_n)


    #values taken from charmm
    AIJ = 3.15
    BIJ = 1.80

    # individual values taken frmo herschbach and laurie 1961
    if (p1 < HELIUM)
    
        if      (p2 < HELIUM) ; AIJ = 1.26; BIJ = 0.025;  # 0.025 is not an error!
        elseif (p2 < NEON)   ; AIJ = 1.66; BIJ = 0.30; 
        elseif (p2 < ARGON)  ; AIJ = 1.84; BIJ = 0.38; 
        elseif (p2 < KRYPTN) ; AIJ = 1.98; BIJ = 0.49; 
        elseif (p2 < XENON)  ; AIJ = 2.03; BIJ = 0.51; 
        elseif (p2 < RADON)  ; AIJ = 2.03; BIJ = 0.25; end

    elseif (p1 < NEON)
        
        if 	    (p2 < NEON)  ;  AIJ = 1.91; BIJ = 0.68; 
        elseif (p2 < ARGON) ;  AIJ = 2.28; BIJ = 0.74; 
        elseif (p2 < KRYPTN);  AIJ = 2.35; BIJ = 0.85; 
        elseif (p2 < XENON) ;  AIJ = 2.33; BIJ = 0.68; 
        elseif (p2 < RADON) ;  AIJ = 2.50; BIJ = 0.97; end
        
    elseif (p1 < ARGON)
        
        if 		(p2 < ARGON) ;  AIJ = 2.41; BIJ = 1.18; 
        elseif (p2 < KRYPTN);  AIJ = 2.52; BIJ = 1.02; 
        elseif (p2 < XENON) ;	AIJ = 2.61; BIJ = 1.28; 
        elseif (p2 < RADON) ;	AIJ = 2.60; BIJ = 0.84; end
        
    elseif (p1 < KRYPTN)
        
        if 		(p2 < KRYPTN) ; AIJ = 2.58; BIJ = 1.41;
        elseif (p2 < XENON)  ; AIJ = 2.66; BIJ = 0.86;
        elseif (p2 < RADON)  ; AIJ = 2.75; BIJ = 1.14; end
        
    elseif (p1 < XENON)
        
        if 		(p2 < XENON) ; AIJ = 2.85; BIJ = 1.62;
        elseif (p2 < XENON) ; AIJ = 2.76; BIJ = 1.25; end
        
    end
    kb = ((AIJ - BIJ) / (r0 - BIJ))^3
    return kb
end 


function setupStretches(msb::MStretchBendComponent{T}, aro_rings) where T<:Real
    stretch_params     = msb.ff.parameters.sections["BondStretch"].data
    stretch_e_params   = msb.ff.parameters.sections["EmpiricalBondParams"].data
    stretches = Dict{NTuple{2, Atom{T}} , MStretch}()
    gdf_stretch = groupby(stretch_params, [:b, :c])
    gdf_stretch_e = groupby(stretch_e_params, [:type1, :type2])
    bonds = non_hydrogen_bonds(msb.ff.system)
    unassigned_atoms = Atom{T}[]

    for bond in bonds
        at1 = atom_by_idx(msb.ff.system, bond.a1)
        at2 = atom_by_idx(msb.ff.system, bond.a2)
        #at1.idx > at2.idx && begin at1, at2 = at2, at1 end
        at1_t = parse(Int8, at1.atom_type)
        at2_t = parse(Int8, at2.atom_type)
        is_sbmb = has_property(bond, :MMFF94SBMB)
        row = get(gdf_stretch, (b = at1_t, c = at2_t), missing)
        
        if !ismissing(row)
            
            stretches[Tuple(sort([at1, at2], by=at->at.idx))] = MStretch(at1, at2, 
                T(row.kb[1 + is_sbmb]), T(row.r0[1 + is_sbmb]), zero(T), is_sbmb, false)
            
            set_property!(bond, :MMFF94RBL, row.r0[1 + is_sbmb])
            continue
        end

        r0 = calculate_stretch_r0(
            bond,
            at1,
            at2,
            msb.ff.parameters.radii,
            msb.ff.parameters.electronegativities,
            msb.ff.parameters.sections["AtomTypeProperties"].data, 
            aro_rings     
        )

        kb = calculate_stretch_kb(r0, at1, at2, gdf_stretch_e)

        if r0 != -1 && kb != -1
            stretches[Tuple(sort([at1, at2], by=at->at.idx))] = MStretch(at1, at2,
                T(kb), T(r0), zero(T), is_sbmb, true)
            set_property!(bond, :MMFF94RBL, T(r0))
            continue
        end
 
        @info """Cannot find stretch parameters for atom types $(at1.atom_type) $(at2.atom_type).
        \nAtoms are: \n$((at1.idx, at1.element, at1.atom_type, at1.name)),
        \n$((at2.idx, at2.element, at2.atom_type, at2.name))"""
        push!(unassigned_atoms, at1, at2)
        
    end

    msb.stretches = stretches
    return unassigned_atoms::Vector{Atom{T}}
end

function setupBends(msb::MStretchBendComponent{T}, all_rings) where T<:Real
    atom_type_properties = msb.ff.parameters.sections["AtomTypeProperties"].data
    bend_params          = msb.ff.parameters.sections["AngleBend"].data
    equivs               = msb.ff.parameters.sections["Equivalences"].data
    gdf_b = groupby(bend_params, [:a, :b, :c, :d])
    flat_rings = collect(Iterators.flatten(all_rings))

    unassigned_atoms = Atom{T}[]
    for atom in atoms(msb.ff.system)
        for b1 in non_hydrogen_bonds(atom)
            for b2 in non_hydrogen_bonds(atom)[2:end]
                b1.idx == b2.idx && continue
                at1 = get_partner(b1, atom)
                at2 = atom
                at3 = get_partner(b2, atom)

                if parse(Int8,at1.atom_type) > parse(Int8,at3.atom_type)
                    at1, at3 = at3, at1
                end

                ATIJK = get_bend_type(b1, b2, at1, at2, at3, flat_rings)
                is_lin = Bool(atom_type_properties[at_idx(at2.atom_type), :lin])


                ka, theta0 = get_bend_params(ATIJK, at1.atom_type, at2.atom_type,
                    at3.atom_type, gdf_b, equivs)

                if !ismissing(ka)
                    #l235
                    if ka == 0
                        ka = calculate_bend_empirical_force_constant(
                            at1, at2, at3, theta0, flat_rings
                        )
                    end

                    if ka > 0
                        push!(msb.bends, MBend{T}( theta0, 0, 0, ka, at1, at2, at3,
                            is_lin, ATIJK, Vector3(0,0,0), Vector3(0,0,0)))
                        continue
                    end
                end

                #l250
                theta0 = calculate_bend_empirical_reference_angle(at1, at2, at3, bend_params, flat_rings)
                ka = calculate_bend_empirical_force_constant(at1, at2, at3, theta0, flat_rings)


                if ka != -1 
                    bend = MBend{T}(theta0, 0., 0., ka, at1, at2, at3, is_lin, ATIJK,
                    Vector3(0,0,0), Vector3(0,0,0))
                    push!(msb.bends, bend)
                    continue
                end

                push!(unassigned_atoms, at1, at2, at3)

                @info """MStretchBendComponent - setupBends: Cannot find bend parameters
                for atom types a1:$(at1.atom_type) a2:$(at2.atom_type) a3:$(at1.atom_type)
                    and bend type $ATIJK.\n Atoms are:\n
                    $((at1.idx, at1.element, at1.atom_type, at1.name)), 
                    $((at2.idx, at2.element, at2.atom_type, at2.name)), 
                    $((at3.idx, at3.element, at3.atom_type, at3.name))
                    """
            end
        end
    end

    return
end

function calculate_bend_empirical_reference_angle(at1, at2, at3, atom_type_properties, flat_rings)
    ring = findfirst(
        ring -> 3 <= length(ring) <= 4 && all(at -> at in ring, (at1, at2, at3)),
        flat_rings
        )

    !isnothing(ring) && length(flat_rings[ring]) == 3 ? (return 60) : (return 90)
    atom_props = atom_type_properties[at_idx(at1.atom_type), :]
    # check atom_props fieldnames

    atom_props.crd == 4 && return 109.45
    if      atom_props.crd == 2
        at2.element == Elements.O && return 105.
        Int(at2.element) > 10 && return 95.
        Bool(atom_props.lin) && return 180.
    elseif  atom_props.crd == 3 && atom_props.val == 3 && atom_props.mltb == 0
        at2.element == Elements.N && return 107.
        return 92.
    end

    #default value
    return 120.
end



function calculate_bend_empirical_force_constant(at1, at2, at3, theta0, flat_rings)
    #deg_to_rad
    deg_to_rad = π/180
    theta0 *= deg_to_rad
    at_1_idx = findfirst(==(at1.element), bend_elems)
    at_3_idx = findfirst(==(at2.element), bend_elems)
    at_2_idx = findfirst(==(at3.element), bend_elems)
    any(isnothing, (at_1_idx, at_2_idx, at_3_idx)) && return -1.

    zcz = bend_z_[at_1_idx] * bend_c_[at_2_idx] * bend_z_[at_3_idx]
    zcz == 0 && return -1.

    beta = 1.75
    #check if all 3 atoms in a ring of size3 or 4
    ring = findfirst(
            ring -> 3 <= length(ring) <= 4 && all(at -> at in ring, (at1, at2, at3)), 
            flat_rings
        )

    !isnothing(ring) && length(flat_rings[ring]) == 4 && (beta *= 0.85)
    !isnothing(ring) && length(flat_rings[ring]) == 3 && (beta *= 0.05)

    r01_b = findfirst(b -> b.a1 in (at1.idx, at2.idx) && b.a2 in (at1.idx, at2.idx), bonds(at1))
    r02_b = findfirst(b -> b.a1 in (at2.idx, at3.idx) && b.a2 in (at2.idx, at3.idx), bonds(at2))
    r01 = get_property(bonds(at1)[r01_b], :MMFF94RBL)
    r02 = get_property(bonds(at2)[r02_b], :MMFF94RBL)

    rr = r01 + r02
    D = (r01 - r02)^2 / (r01 + r02)^2
    ex = exp(-2*D)
    asq = theta0^(-2)
    k = beta * zcz * asq * ex / (rr)
    return k

end

function get_bend_params(b_t, a1_t, a2_t, a3_t, gdf, equivs)
    row = coalesce(
        get(gdf, (a=b_t, b=       at_idx(a1_t)     , c=at_idx(a2_t), d=       at_idx(a3_t)     ), missing),
        get(gdf, (a=b_t, b=equivs[at_idx(a1_t), :b], c=at_idx(a2_t), d=equivs[at_idx(a3_t), :b]), missing),
        get(gdf, (a=b_t, b=equivs[at_idx(a1_t), :c], c=at_idx(a2_t), d=equivs[at_idx(a3_t), :c]), missing),
        get(gdf, (a=b_t, b=equivs[at_idx(a1_t), :d], c=at_idx(a2_t), d=equivs[at_idx(a3_t), :d]), missing),
        get(gdf, (a=b_t, b=equivs[at_idx(a1_t), :e], c=at_idx(a2_t), d=equivs[at_idx(a3_t), :e]), missing),
    )

    if !ismissing(row)
        return only(row.ka), only(row.theta0)
    end

    return missing, missing

end

function get_bend_type(b1, b2, at1, at2, at3, flat_rings)
    # this function assumes that an atom can't be in a ring of three members
    # and a ring of 4 members t the same time
    ring = findfirst(
            ring -> 3 <= length(ring) <= 4 && all(at -> at in ring, (at1, at2, at3)),
            flat_rings
        )
        

    sum_bond_types = sum(has_property(b1, :MMFF94SBMB) + has_property(b2, :MMFF94SBMB))
    isnothing(ring) && return sum_bond_types

    if length(flat_rings[ring]) == 3
        sum_bond_types != 0 ? (return 4 + sum_bond_types) : (return 3)
    end

    if length(flat_rings[ring]) == 4
        sum_bond_types != 0 ? (return 6 + sum_bond_types) : (return 4)
    end

end

function setupStretchBends(msb::MStretchBendComponent{T}) where T<:Real
    gdf = groupby(msb.ff.parameters.sections["StretchBend"].data, [:type, :I, :J, :K])
    gdf2 = groupby(msb.ff.parameters.sections["StretchBendEmpirical"].data, [:IR, :JR, :KR])

    for bend in msb.bends

        stretch1 = msb.stretches[Tuple(sort([bend.at1, bend.at2], by=at -> at.idx))]
        stretch2 = msb.stretches[Tuple(sort([bend.at2, bend.at3], by=at -> at.idx))]
        
        SBTIJK = calculate_SBTIJK(bend.ATIJK, stretch1.sbmb)
        SBTIJK == 2 && bend.at1.atom_type == bend.at2.atom_type && (SBTIJK = 1)

        kba_ijk, kba_kji = get_kba(gdf, gdf2, SBTIJK, bend.at1, bend.at2, bend.at3)
        ismissing(kba_ijk) && @info "Error in StretchBend setup for $(bend.at1.idx), $(bend.at2.idx), $(bend.at3.idx)"
        ismissing(kba_ijk) && continue

        push!(msb.stretch_bends, MStretchBend{T}(kba_ijk, kba_kji, SBTIJK,
            stretch1, stretch2, bend))
    end 

end

function get_period(at)
    Int8(at.element) in 1:2   && return 0
    Int8(at.element) in 3:10  && return 1
    Int8(at.element) in 11:18 && return 2
    Int8(at.element) in 19:36 && return 3
    Int8(at.element) in 37:54 && return 4
    return 0
end

function get_kba(gdf, gdf2, SBTIJK, at1, at2, at3)
    row = get(gdf,
        (type=SBTIJK, I=parse(Int8,at1.atom_type), J=parse(Int8,at2.atom_type),
            K=parse(Int8,at3.atom_type)),
        missing
    )

    !ismissing(row) && return only(row)[[:kbaIJK, :kbaKJI]]

    #is ismissing then index by atom period
    #@info "was missing; $(at1.atom_type) $(at2.atom_type) $(at3.atom_type) using elements $(at1.idx), $(at2.idx), $(at3.idx)"

    
    row = coalesce(
        get(gdf2, (IR=get_period(at1), JR=get_period(at2), KR=get_period(at3)), missing),
        get(gdf2, (IR=get_period(at3), JR=get_period(at2), KR=get_period(at1)), missing)
    )


    if !ismissing(row) 
        get_period(at1) <= get_period(at3) ? 
            (return only(row)[[:fijk, :kji]]) :
            (return only(row)[[:kji, :fijk]])
    end

    x = 1
    return missing, missing
end


#= "The stretch-bend types are defined in terms of the constituent bond types BTIJ 
		 and BTJK and the angle type ATIJK as shown below:"

	   Stretch-Bend       Angle           ---- Bond Type ----
					 Type         Type            I-J             J-K
	 -----------------------------------------------------------
						0             0               0               0
						1             1               1               0
						2             1(*)            0               1 		* error in CHARMM docu
						3             2               1               1
						4             4               0               0
						5             3               0               0
						6             5               1               0
						7             5               0               1
						8             6               1               1
						9             7               1               0
					 10             7               0               1
					 11             8 (*)           1               1
=#
function calculate_SBTIJK(angle_type, b_type1)
    angle_type == 0 && return angle_type
    angle_type == 1 && b_type1 ? (return 1) : (return 2)
    angle_type == 2 && return 3
    angle_type == 3 && return 5
    angle_type == 4 && return angle_type
    angle_type == 5 && b_type1 ? (return 6) : (return 7)
    angle_type == 6 && return 8
    angle_type == 7 && b_type1 ? (return 9) : (return 10)
    angle_type == 8 && return 11

    @info "Could not calculate SBTIJK for $angle_type, $b_type1, $b_type2"

    return -1

end