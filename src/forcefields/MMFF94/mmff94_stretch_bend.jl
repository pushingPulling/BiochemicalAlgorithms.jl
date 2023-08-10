@auto_hash_equals mutable struct MBend{T<:Real}
    theta0   ::T            #in degrees
    theta_d  ::T            #in degrees
    theta    ::T            #in degrees
    ka       ::T            #in millidyne * angstrom / rad^2
    at1      ::Atom{T}
    at2      ::Atom{T}
    at3      ::Atom{T}
    is_linear::Bool
    ATIJK    ::Int
    n1       ::Vector3{T}
    n2       ::Vector3{T}
end

@auto_hash_equals mutable struct MStretch{T<:Real}
    at1      ::Atom{T}
    at2      ::Atom{T}
    kb       ::T            #in millidyne / angstrom
    r0       ::T            #in angstrom
    delta_r  ::T            #in angstrom
    sbmb     ::Bool
    empirical::Bool
    n        ::Vector3{T}
end

@auto_hash_equals mutable struct MStretchBend{T<:Real}
    kba_ijk     ::T         #in millidyne /rad
    kba_kji     ::T         #in millidyne /rad
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
    rings_cache     ::Tuple{Vector{Vector{Atom{T}}}, Vector{Vector{Atom{T}}}}
end

function MStretchBendComponent{T}(ff, all_rings, aromatic_rings) where T<:Real
    all_rings_ = collect(Iterators.flatten(all_rings))
    aromatic_rings_ = collect(Iterators.flatten(aromatic_rings))

    comp = MStretchBendComponent{T}(ff, "msbc", Dict{Symbol, Any}(),
        Dict{String, T}(), Dict{NTuple{2,Atom{T}}, MStretch{T}}(), 
        MBend{T}[], MStretchBend{T}[], (all_rings_, aromatic_rings_))

    setup!(comp)
    update!(comp)
    comp
end


function update!(msb::MStretchBendComponent{T}) where T<:Real
    nothing
end


function setup!(msb::MStretchBendComponent{T}) where T<:Real
    #add ring calc here
    empty!(msb.cache)
    empty!(msb.energy)
    empty!(msb.stretches)
    empty!(msb.bends)
    empty!(msb.stretch_bends)

    unassigned_atoms   = Atom{T}[]
    push!(unassigned_atoms, setupStretches(msb)...)
    setupBends(msb)
    setupStretchBends(msb)
end

function calculate_stretch_r0(bond, at1, at2, radii,
            electro_neg, atom_type_properties, aro_rings)
    

    a1_n = Int(at1.element)
    a2_n = Int(at2.element)

    # default is same value as in mmff94_parameters.jl::assign_to
    (at1.atom_type == "0" || at2.atom_type == "0") && return -1.0
    (a1_n > 53 || a2_n > 53 || a1_n == 0 || a2_n == 0) && return -1.0

    r1 = radii[a1_n]
    r2 = radii[a2_n]
    r1 == 0.0 || r2 == 0.0 && return -1.0

    if a1_n != 0 && a2_n != 0
        bond.order in (BondOrder.Unknown, BondOrder.Quadruple) && return -1.0
        bo = bond.order
        b1 = atom_type_properties[at_idx(at1.atom_type), :mltb]
        b2 = atom_type_properties[at_idx(at2.atom_type), :mltb]

        if     b1 == b2 == 1;    bo = BondOrder.Quadruple;
        elseif b1  + b2 == 3;    bo = BondOrder.Aromatic;
        else
            if any(r -> at1 in r && at2 in r, aro_rings)
                !Bool(atom_type_properties[at_idx(at1.atom_type), :pilp]) &&
                    !Bool(atom_type_properties[at_idx(at2.atom_type), :pilp]) ?
                        bo = BondOrder.Quadruple : bo = BondOrder.Aromatic
            end
        end

        if     bo == BondOrder.Aromatic; r1 -= 0.04 ;    r2 -= 0.04 ;
        elseif bo == BondOrder.Quadruple; r1 -= 0.075;    r2 -= 0.075;
        elseif bo == BondOrder.Triple; r1 -= 0.17 ;    r2 -= 0.17 ;
        elseif bo == BondOrder.Double; r1 -= 0.1  ;    r2 -= 0.1  ;
        else  #bond.order == Single

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
    r0 = r1 + r2 - c * delta_e^n - d
    return r0 

end


function calculate_stretch_kb(r0, at1, at2, gdf_stretch_e)
    a1_n, a2_n = extrema((Int(at1.element), Int(at2.element)))
    row = get(gdf_stretch_e, (type1=a1_n, type2=a2_n), missing)
    if !ismissing(row) 
        kb = only(row.kb) * (only(row.r0)/r0)^6 
        return kb
    end

    at1_n = Int8(at1.element)
    at2_n = Int8(at2.element)
    p1, p2 = extrema((at1_n, at2_n))


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
        # this is consistent with the C++ code. Should it maybe be p2 > XENON instead?
        elseif (p2 < XENON) ; AIJ = 2.76; BIJ = 1.25; end
        
    end
    kb = ((AIJ - BIJ) / (r0 - BIJ))^3
    return kb
end 


function setupStretches(msb::MStretchBendComponent{T}) where T<:Real
    aro_rings          = msb.rings_cache[2]
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
            row_ = row[min(1 + is_sbmb, size(row)[1]), :]
            stretches[Tuple(sort([at1, at2], by=at->at.idx))] = MStretch(at1, at2, 
                T(row_.kb), T(row_.r0), zero(T), is_sbmb, false,
                Vector3{T}(0,0,0))
            
            set_property!(bond, :MMFF94RBL, row_.r0)
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
                T(kb), T(r0), zero(T), is_sbmb, true, Vector3{T}(0,0,0))
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

function setupBends(msb::MStretchBendComponent{T}) where T<:Real
    all_rings            = msb.rings_cache[1]
    atom_type_properties = msb.ff.parameters.sections["AtomTypeProperties"].data
    bend_params          = msb.ff.parameters.sections["AngleBend"].data
    equivs               = msb.ff.parameters.sections["Equivalences"].data
    gdf_b = groupby(bend_params, [:a, :b, :c, :d])

    unassigned_atoms = Atom{T}[]
    for atom in atoms(msb.ff.system)
        for (cur_bond, b1) in enumerate(non_hydrogen_bonds(atom))
            for b2 in non_hydrogen_bonds(atom)[cur_bond + 1:end]
                #b1.idx == b2.idx && continue
                at1 = get_partner(b1, atom)
                at2 = atom
                at3 = get_partner(b2, atom)

                if parse(Int8,at1.atom_type) > parse(Int8,at3.atom_type)
                    at1, at3 = at3, at1
                end

                ATIJK = get_bend_type(b1, b2, at1, at2, at3, all_rings)
                #ToDo fix this
                is_lin = at2.atom_type != "0" ? Bool(atom_type_properties[at_idx(at2.atom_type), :lin]) : false



                ka, theta0 = get_bend_params(ATIJK, at1.atom_type, at2.atom_type,
                    at3.atom_type, gdf_b, equivs)

                if !ismissing(ka)
                    #l235
                    if ka == 0
                        ka = calculate_bend_empirical_force_constant(
                            at1, at2, at3, theta0, all_rings
                        )
                    end

                    if ka > 0
                        push!(msb.bends, MBend{T}( theta0, 0, 0, ka, at1, at2, at3,
                            is_lin, ATIJK, Vector3(0,0,0), Vector3(0,0,0)))
                        continue
                    end
                end

                # expressed in degrees
                theta0 = calculate_bend_empirical_reference_angle(at1, at2,
                    at3, bend_params, all_rings)
                ka = calculate_bend_empirical_force_constant(at1, at2, at3, theta0, all_rings)


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

# all angles in degrees
function calculate_bend_empirical_reference_angle(at1, at2, at3, atom_type_properties, all_rings)
    ring = findfirst(
        ring -> 3 <= length(ring) <= 4 && all(at -> at in ring, (at1, at2, at3)),
        all_rings
    )

    !isnothing(ring) && length(all_rings[ring]) == 3 ? (return 60) : (return 90)
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



function calculate_bend_empirical_force_constant(at1, at2, at3, theta0, all_rings)
    #deg_to_rad
    theta0 = theta0 *1u"°" |> u"rad" |> ustrip
    at_1_idx = findfirst(==(at1.element), bend_elems)
    at_2_idx = findfirst(==(at2.element), bend_elems)
    at_3_idx = findfirst(==(at3.element), bend_elems)
    any(isnothing, (at_1_idx, at_2_idx, at_3_idx)) && return -1.

    r01_b = findfirst(b -> b.a1 in (at1.idx, at2.idx) && b.a2 in (at1.idx, at2.idx), bonds(at1))
    r02_b = findfirst(b -> b.a1 in (at2.idx, at3.idx) && b.a2 in (at2.idx, at3.idx), bonds(at2))
    r01 = get(bonds(at1)[r01_b].properties, :MMFF94RBL, missing)
    r02 = get(bonds(at2)[r02_b].properties, :MMFF94RBL, missing)
    (ismissing(r01) || ismissing(r02)) && return -1

    zcz = bend_z_[at_1_idx] * bend_c_[at_2_idx] * bend_z_[at_3_idx]
    zcz == 0 && return -1.

    beta = 1.75
    #check if all 3 atoms in a ring of size3 or 4
    ring = findfirst(
            ring -> 3 <= length(ring) <= 4 && all(at -> at in ring, (at1, at2, at3)), 
            all_rings
        )

    !isnothing(ring) && length(all_rings[ring]) == 4 && (beta *= 0.85)
    !isnothing(ring) && length(all_rings[ring]) == 3 && (beta *= 0.05)

    r01_b = findfirst(b -> b.a1 in (at1.idx, at2.idx) && b.a2 in (at1.idx, at2.idx), bonds(at1))
    r02_b = findfirst(b -> b.a1 in (at2.idx, at3.idx) && b.a2 in (at2.idx, at3.idx), bonds(at2))

    rr = r01 + r02
    D = (r01 - r02)^2 / (r01 + r02)^2
    ex = exp(-2*D)
    asq = theta0^(-2)
    k = beta * zcz * asq * ex / (rr)

    return k
end

function get_bend_params(b_t, a1_t, a2_t, a3_t, gdf, equivs)
    
    #first row may contain zeros in column :b or :d
    row = get(gdf, (a=b_t, b=       at_idx(a1_t)     , c=at_idx(a2_t), d=       at_idx(a3_t)     ), missing)
    !ismissing(row) && return only(row.ka), only(row.theta0)
    any(==("0"), (a1_t, a3_t)) && return missing, missing
     

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

function get_bend_type(b1, b2, at1, at2, at3, all_rings)
    # this function assumes that an atom can't be in a ring of three members
    # and a ring of 4 members t the same time
    ring = findfirst(
        ring -> 3 <= length(ring) <= 4 && all(at -> at in ring, (at1, at2, at3)),
        all_rings
    )
        

    sum_bond_types = has_property(b1, :MMFF94SBMB) + has_property(b2, :MMFF94SBMB)
    isnothing(ring) && return sum_bond_types

    if length(all_rings[ring]) == 3
        sum_bond_types != 0 ? (return 4 + sum_bond_types) : (return 3)
    end

    if length(all_rings[ring]) == 4
        sum_bond_types != 0 ? (return 6 + sum_bond_types) : (return 4)
    end

end

function setupStretchBends(msb::MStretchBendComponent{T}) where T<:Real
    gdf = groupby(msb.ff.parameters.sections["StretchBend"].data, [:type, :I, :J, :K])
    gdf2 = groupby(msb.ff.parameters.sections["StretchBendEmpirical"].data, [:IR, :JR, :KR])

    for bend in (b for b in msb.bends if !b.is_linear)

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


function calculate_delta(bend::MBend{T}) where T<:Real
    v1 = bend.at1.r - bend.at2.r
    v2 = bend.at3.r - bend.at2.r

    v1_len = norm(v1)
    v2_len = norm(v2)

    v1 /= v1_len
    v2 /= v2_len

    if any(isnan, first.((v1_len, v2_len)))
        bend.n1 = Vector3{T}(0., 0., 0.)
        bend.n2 = Vector3{T}(0., 0., 0.)
        return 
    end

    c = normalize(cross(v1, v2))

    if isnan(first(c))
        bend.n1 = Vector3{T}(0., 0., 0.)
        bend.n2 = Vector3{T}(0., 0., 0.)
        return 
    end

    costheta = dot(v1, v2)
    if      costheta >  1.0;      theta = 0.        ; costheta =  1.;
    elseif  costheta < -1.0;      theta = Unitful.pi; costheta = -1.;
    else                          theta = acos(costheta)
    end

    if bend.is_linear
        bend.n1 = (v2 - (v1 * costheta)) / v1_len
        bend.n2 = (v1 + (v2 * costheta)) / v2_len
    end

    theta = theta * 1u"rad" |> u"°" |> ustrip

    bend.theta_d = theta - bend.theta0
    bend.theta = theta

    if !bend.is_linear
        t1 = normalize(cross(v1, c))
        t2 = normalize(cross(v2, c))

        bend.n1 = -t1 / v1_len
        bend.n2 =  t2 / v2_len
    end

end


function calculate_delta(stretch::MStretch{T}) where T<:Real
    direction = stretch.at1.r - stretch.at2.r
    distance = norm(direction)
    stretch.delta_r = distance - Float64(stretch.r0)
    stretch.n = direction / distance
end


function compute_forces(msb::MStretchBendComponent{T}) where T<:Real

    map(calculate_delta, msb.bends)
    map(calculate_delta, values(msb.stretches))

    get(msb.ff.options, :bends_enabled, true)           && map(compute_bend_forces, msb.bends)
    get(msb.ff.options, :stretches_enabled, true)        && map(compute_stretch_forces, values(msb.stretches))
    get(msb.ff.options, :stretch_bends_enabled, true)  && map(compute_stretch_bend_forces, msb.stretch_bends)

    nothing
end

function compute_energy(msb::MStretchBendComponent{T}) where T<:Real

    map(calculate_delta, msb.bends)
    map(calculate_delta, values(msb.stretches))

    energy = 0

    get(msb.ff.options, :bends_enabled, true)          && (energy += mapreduce(compute_bend_energy, +, msb.bends,init=zero(T)))
    get(msb.ff.options, :stretches_enabled, true)      && (energy += mapreduce(compute_stretch_energy, +, values(msb.stretches),init=zero(T)))
    get(msb.ff.options, :stretch_bends_enabled, true)  && (energy += mapreduce(compute_stretch_bend_energy, +, msb.stretch_bends,init=zero(T)))
    return energy
end

function compute_bend_forces(bend::MBend{T}) where T<:Real
    factor = bend.is_linear ? -bend_kx * bend.ka : 
        bend_k₀_forces * bend.ka * bend.theta_d * (2. + 3. * bend_k1 * bend.theta_d)

    n1 = bend.n1 * factor * (1000 * 1e10/Constants.N_A_nounit)
    n2 = bend.n2 * factor * (1000 * 1e10/Constants.N_A_nounit)

    bend.at1.F += n1
    bend.at3.F += n2

    bend.at2.F -= n1
    bend.at2.F -= n2
end

function compute_stretch_forces(stretch::MStretch{T}) where T<:Real
    a = stretch_k0 * stretch.kb * stretch.delta_r
    force = (2 + 3 * stretch_cubic_strength_constant * stretch.delta_r + 
        4 * stretch_kcs * stretch.delta_r^2) * a
    
    direction = stretch.n * ((force)u"kJ/mol/Å"/Constants.N_A |> u"N" |> ustrip)
    stretch.at1.F -= direction
    stretch.at2.F += direction
end


function compute_stretch_bend_forces(sb::MStretchBend{T}) where T<:Real
    bend = sb.bend
    c1 = normalize(bend.at1.r - bend.at2.r)
    c2 = normalize(bend.at3.r - bend.at2.r)

    temp = bend.theta_d * stretch_bend_k0 * force_prefactor
    r1 = c1 * sb.kba_ijk * temp
    r3 = c2 * sb.kba_kji * temp

    d_ij = sb.stretch_i_j.delta_r
    d_jk = sb.stretch_j_k.delta_r

    sb_scale = -(sb.kba_ijk * d_ij + sb.kba_kji * d_jk) * sb_constant
    r1 += bend.n1 * sb_scale
    r3 += bend.n2 * sb_scale

    bend.at1.F -= r1
    bend.at2.F += r1 + r3
    bend.at3.F -= r3
end


function compute_stretch_bend_energy(sb::MStretchBend{T}) where T<:Real
    bend = sb.bend
    
    sb1 = stretch_bend_k0 * (sb.kba_ijk * sb.stretch_i_j.delta_r) * bend.theta_d
    sb2 = stretch_bend_k0 * (sb.kba_kji * sb.stretch_j_k.delta_r) * bend.theta_d
    sb1 + sb2
end

function compute_stretch_energy(stretch::MStretch{T}) where T<:Real
    delta = stretch.delta_r

    eb_ij = stretch_k0 * stretch.kb * delta^2 * 
        (Float64(1) + 
            stretch_cubic_strength_constant * delta + 
            Float64(stretch_kcs) * delta^2)

end

function compute_bend_energy(bend::MBend{T}) where T<:Real
    energy = (bend.is_linear 
        ? bend_kx * bend.ka * (1. + cos(deg2rad(bend.theta))) 
        : bend_k₀ * bend.ka * bend.theta_d^2 * (1. + bend_k1 * bend.theta_d))

    #bend.energy = energy
end
