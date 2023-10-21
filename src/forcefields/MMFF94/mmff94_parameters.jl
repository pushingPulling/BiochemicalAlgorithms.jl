using AutoHashEquals
using CSV
using Pipe
export MMFF94Parameters, _read_types_paramfile, assign_atom_types_mmff94,
assign_bond_types_mmff94, assign_charges

###### Constants #######
bend_elems = (Elements.H, Elements.B,Elements.C,Elements.N,Elements.O,Elements.F,
    Elements.Si,Elements.P,Elements.S,Elements.Cl,Elements.Br,Elements.As,Elements.I)
bend_z_ = (1.395, 0. , 2.494, 2.711, 3.045, 2.847, 2.350, 2.350, 2.980, 2.909, 3.017, 0., 3.086)
bend_c_ = (0.,    0.704, 1.016, 1.113, 1.337, 0., 0.811, 1.068, 1.249, 1.078, 0., 0.825, 0.)

HELIUM    = 2
NEON      = 10
ARGON     = 18
KRYPTN    = 36
XENON     = 54
RADON     = 86
########################

#helperfunc to index AtomTypeProperties in forcefield parameters
at_idx(type) = parse(Int8, type) > 82 ? parse(Int8, type) - 4 : parse(Int8, type)


@auto_hash_equals struct MMFF94Parameters{T<:Real} <: AbstractForceFieldParameters
    sections::OrderedDict{String, BALLIniFileSection}
    radii::Vector{T}
    electronegativities::Vector{T}

    function MMFF94Parameters{T}(
            path::String = ball_data_path("forcefields/MMFF94/mmff94.ini")) where {T<:Real}

        ini_file = read_ball_ini_file(path, T)

        sections = ini_file.sections
        
        # see http://www.ccl.net/cca/data/MMFF94/
        radii = T.([
             0.33, 0.0,
             1.34, 0.90, 0.81, 0.77, 0.73, 0.72, 0.74, 0.0,
             1.54, 1.30, 1.22, 1.15, 1.09, 1.03, 1.01, 0.0,
             1.96, 1.74,
             1.44, 1.36, 0.00, 0.00, 0.00,
             0.00, 0.00, 0.00, 1.38, 1.31,
             1.19, 1.20, 1.20, 1.16, 1.15, 0.0,
             2.11, 1.92,
             1.62, 1.48, 0.00, 0.00, 0.00,
             0.00, 0.00, 0.00, 1.53, 1.48,
             1.46, 1.40, 1.41, 1.35, 1.33, 0.0
        ])
        
        electronegativities = T.([
             2.20, 0.0,
             0.97, 1.47, 2.01, 2.5, 3.07, 3.5, 4.10, 0.0,
             1.01, 1.23, 1.47, 1.74, 2.06, 2.44, 2.83, 0.0,
             0.91, 1.04,
             1.3, 1.5, 1.6, 1.6, 1.5,
             1.8, 1.8, 1.8, 1.9, 1.6,
             1.82, 2.02, 2.20, 2.48, 2.74, 0.0,
             0.89, 0.99,
             1.3, 1.4, 1.6, 1.8, 1.9,
             2.2, 2.2, 2.2, 1.9, 1.7,
             1.49, 1.72, 1.82, 2.01, 2.21, 0.0
        ])

        new(sections, radii, electronegativities)
    end
end

function MMFF94Parameters(path::String = ball_data_path("forcefields/MMFF94/mmff94.ini"))
    MMFF94Parameters{Float32}(path)
end


Base.show(io::IO, mmff_param::MMFF94Parameters{T}) where {T<:Real} = 
    print(io, 
        "MMFF94 parameters with $(length(mmff_param.sections)) sections.")


struct AtomTypeParams
    symbol::String7
    type::String3
    rule::String127
end
        
Base.show(io::IO, x::AtomTypeParams) = show(io::IO, (x.symbol, x.type, x.rule))

function _read_types_paramfile(path)
    lines = filter!(line -> !startswith(line, "*") ,readlines(path))
    lines = split.(lines, r"\s+\|\s+")
    types_param_dict = Dict{String3, Vector{AtomTypeParams}}()
    for line in lines
        vec = get!(types_param_dict, String3(line[1]), Vector{AtomTypeParams}())  
        push!(vec, AtomTypeParams(line[2:4]...))
    end
    return types_param_dict
end
        
        
        
function assign_atom_types_mmff94(ac::AbstractAtomContainer, mmff94_params, aro_rings)
    #add aro rings, add atomtype properties
    #read contents of TYPES.PAR
    path = ball_data_path("forcefields/MMFF94/TYPES.PAR")
    elem_to_params = _read_types_paramfile(path)

    ######################
    # setup aromatic
    ######################

    df = mmff94_params.sections["AtomTypeProperties"].data
    hetero_atom_types_df = @rsubset df begin
        :type != 10 && 
        (
            (:aspec == 7 && :crd == :val == 3) ||
            ((:aspec == 8 || :aspec == 16) && :crd == :val == 2)
        )
        @kwarg view = true
    end
    
    df = mmff94_params.sections["Aromatic"].data
    cation_df = @rsubset df begin
        !(endswith(:oldtype, "*"))
        :rsize == 5
        Bool(:imcat)
        @kwarg view = true
    end

    ring_5_df = @rsubset df begin
        !(endswith(:oldtype, "*"))
        :rsize == 5
        @kwarg view = true
    end
    
    
    #map(mol -> assign_to(mol, elem_to_params), molecules(ac))
    for mol in molecules(ac)
        assign_to(mol, elem_to_params)
    end

    ##############
    # assign heterocyclic 5 ring member types
    # e.g. CSA CSB N5A N5B C5
    ##############


    assign_heterocyclic_member_types(aro_rings[1],
        cation_df, 
        ring_5_df, 
        mmff94_params.sections["Symbols"].data,
        hetero_atom_types_df
    )


    ## assign hydrogens
    assign_hydrogen_types(ac,
        mmff94_params.sections["HydrogenTypes"].data,
        mmff94_params.sections["Symbols"].data
    )

    untyped_atoms_count = count(name -> name == "Any", atoms_df(ac).name)
    if untyped_atoms_count > 0
        @info "Could not assign atom type for $untyped_atoms_count atoms."
    end
end


function assign_hydrogen_types(ac, h_df, sym_df)
    hydro_adjacent_df = filter(at -> at.element == Elements.H && nbonds(at) == 1, atoms(ac))
    for at in hydro_adjacent_df
        partner = get_partner(first(bonds(at)), at)
        hydro_df_idx = findfirst(==(partner.name), h_df.a)
        isnothing(hydro_df_idx) && return
        h_name = h_df.b[hydro_df_idx]
        H_nr = sym_df.type[findfirst(==(h_name), sym_df.symbol)]
        parse(Int8, at.atom_type) > H_nr && continue
        at.name = h_name
        at.atom_type = string(H_nr)
    end
end

function assign_heterocyclic_member_types(aro_rings, cation_df, ring_5_df, symbols_df, hetero_atom_types_df)
    for ring in filter(r -> length(r) == 6, aro_rings)
        for atom in ring
            if atom.name == "CNN+"
                atom.name = "CB"
                atom.atom_type = "37"
                map(bonds(atom)) do bond
                    p = get_partner(bond, atom)
                    p.atom_type == "55" && (p.atom_type = "40"; p.name = "NC=N")
                end
            end
        end
    end

    for ring in filter(r -> length(r) == 5, aro_rings)
        hetero_atom::Union{Missing, Atom} = missing
        cation_atom::Union{Missing, Atom} = missing
        anion_atom ::Union{Missing, Atom} = missing

        for at in ring
            if (
                    !in(at.element, (Elements.C, Elements.N)) || 
                    in(parse(Int, at.atom_type), hetero_atom_types_df.type)
                ) &&
                (
                    ismissing(hetero_atom) ||
                    hetero_atom.element == Elements.N
                )
                hetero_atom = at
            end

            at.name in cation_df.oldtype && (cation_atom = at)
            at.name in ("NM", "N5M") && (anion_atom = at)
        end

        charged = ismissing(hetero_atom) && any(ismissing, (cation_atom, anion_atom))
        charged_atom = ismissing(cation_atom) ? anion_atom   : cation_atom
        atom_x       = ismissing(hetero_atom) ? charged_atom : hetero_atom

        ismissing(atom_x) && continue
        L5 = 0

        for at in ring
            if hetero_atom == at        #test this ToDo
                L5 = 1
            elseif charged
                L5 = 4
            elseif is_bound_to(atom_x, at)
                L5 = 2
            else
                L5 = 3
            end
        
            assign_5_ring_type(
                at,
                L5,
                ismissing(anion_atom),
                ismissing(cation_atom),
                ring_5_df,
                symbols_df
            )
        end

    end 
end
        
function assign_5_ring_type(atom, L5, is_anion, is_cation, ring_5_df, symbols_df)
    oldtype = atom.name
    is_anion && oldtype == "N5M" && return true
    #row = ring_5_df[(ring_5_df.oldtype .== oldtype) .& (ring_5_df.l5 .== L5), :]
    idx = findfirst(x-> x.oldtype == oldtype && x.l5 == L5, eachrow(ring_5_df))
    if !isnothing(idx) && is_cation && is_anion
        row = ring_5_df[idx, :]
        found = true
        !Bool(row.imcat)   && (found = false)
        !Bool(row.n5anion) && (found = false)
        if found
            atom.name = row.aromtype
            atom.atom_type = string(symbols_df[findfirst(x -> x == row.arom_type, symbols_df.symbol), :type])
            return true
        end
    end


    new_type = ""
    if atom.element == Elements.C
        L5 == 2 && (new_type = "C5A" )
        L5 == 3 && (new_type = "C5B" )
        L5 == 4 && (new_type = "C5"  )
    elseif atom.element == Elements.N
        L5 == 1 && (new_type = "NPYL")
        L5 == 2 && (new_type = "N5A" )
        L5 == 3 && (new_type = "N5B" )
        L5 == 4 && (new_type = "N5"  )

        oxide = false
        val   = 0

        for bond in bonds(atom)
            !in(bond.order, (BondOrder.Unknown, BondOrder.Aromatic)) && (val += Int(bond.order))
            bond.order == BondOrder.Aromatic && (val += 5)
            partner = get_partner(bond, atom)
            partner.element == Elements.O && nbonds(atom) == 1 && (oxide = true)
        end

        if val == 4
            if oxide
                new_type *= length(new_type) == 3 ? "X" : "OX"
            else
                new_type *= "+"
            end
        end

    elseif atom.element == Elements.O
        L5 == 1 && (new_type = "OFUR")
    elseif atom.element == Elements.S
        L5 == 1 && (new_type = "STHI")
    else
        @info "In \"mmff94_parameters.jl::assign_5_ring_type\" : No valid MMFF94 Type was assigned for $atom"
    end

    atom.name      = new_type
    atom.atom_type = string(symbols_df[findfirst(x -> x == new_type, symbols_df.symbol), :type])

end

# at this point i have a list of all rings and a list of all aromatic rings
# determines the atom types for all atoms in `mol`. ´all_rings` is a vector
# of rings in this molecule and ´aromatic_rings´ is a vector of aromatic
# rings in that molecule 
function assign_atom_type_decision_tree(at::Atom{T}, all_rings, rings_5, rings_6, rings_aro) where T<:Real

    # this function is potentially recursive; check if type already assigned
    at.atom_type != "0" && (return)
    hydrogen_count = sulphur_count = oxygen_count = nitrogen_count = 0
    bonded_elem = nothing
    is_NNN_or_NNC = false

    #check if at in ring5/6
    # check if in same ring
    # get explicit valence - sum  of bond orders to this atom
    is_in_ring_5(at) = any(ring-> at in ring, rings_5)
    is_in_ring_6(at) = any(ring-> at in ring, rings_6)
    is_aromatic(at) = any(ring-> at in ring, rings_aro)
    in_same_ring(at, p_at) = any(mol -> any(ring-> at in ring && p_at in ring, mol),all_rings)
    neighbors(at) = (get_partner(bond, at) for bond in bonds(at))
    #explicit_degree(atom) = count(p_at -> !(p_at.element == Elements.H && nbonds(p_at) == 1), neighbors(atom)) # prev version: length(nonhydro_bonds)
    explicit_degree(at) = length(non_hydrogen_bonds(at))
    bonded_aromatic(at, at2) = is_aromatic(at) && in_same_ring(at, at2)

    #there are no atoms with Aromatic Order at this point of the program; Therec could be some with order unknown tho
    
    explicit_valence(at) = isempty(bonds(at)) ? nothing : sum(Int.(getproperty.(bonds(at),:order)))
    # explicit_valence(atom) = begin 
    #     partners = filter(p_at -> !(p_at.element == Elements.H && nbonds(p_at) == 1), collect(neighbors(atom))) # prev version: length(nonhydro_bonds)
    #     isempty(partners) ? nothing : sum(Int.(getproperty.((bonded(atom, p) for p in partners),:order)))
    # end

    if is_aromatic(at)
        if is_in_ring_5(at)
            #Aromatic 5-ring sulfur with pi lone pair (STHI)
            if at.element == Elements.S; at.atom_type = "44";return end
            #Aromatic 5-ring oxygen with pi lone pair (OFUR)
            if at.element == Elements.O; at.atom_type = "59";return end
            if at.element == Elements.N
                if any(
                        partner_at -> partner_at.element == Elements.O && explicit_degree(partner_at) == 1,
                        (get_partner(bond, at) for bond in bonds(at))
                    )
                    at.atom_type = "82"; return
                    # N-oxide nitrogen in 5-ring alpha position,
                    # N-oxide nitrogen in 5-ring beta position,
                    # N-oxide nitrogen in other 5-ring  position,
                    # (N5AX, N5BX, N5OX)
                end
            end

            aro_flag = false
            alpha_pos = []; alpha_atoms = []
            beta_pos =  []; beta_atoms  = []
            for partner_at in neighbors(at)
                !(bonded_aromatic(partner_at, at) && is_in_ring_5(partner_at)) && continue
                in_same_ring(at, partner_at) && (push!(alpha_pos, partner_at;))

                for p_partner_at in (get_partner(bond, partner_at) for bond in bonds(partner_at))
                    p_partner_at.idx == at.idx && continue
                    ( !(is_aromatic(p_partner_at)) || !is_in_ring_5(p_partner_at) ) && continue
                    aro_flag = true
                    in_same_ring(p_partner_at, partner_at) && (push!(beta_pos, p_partner_at;))
                end
            end

            if aro_flag
                #continue from l 58
                for a_at in alpha_pos
                    if     a_at.element == Elements.S; push!(alpha_atoms, a_at)
                    elseif a_at.element == Elements.O; push!(alpha_atoms, a_at)
                    elseif a_at.element == Elements.N && explicit_degree(a_at) == 3
                        if !any(p_at -> p_at.element == Elements.O && explicit_degree(p_at) == 1)
                            push!(alpha_atoms, a_at)
                        end
                    end
                end

                for b_at in beta_pos
                    if     b_at.element == Elements.S; push!(beta_atoms, b_at)
                    elseif b_at.element == Elements.O; push!(beta_atoms, b_at)
                    elseif b_at.element == Elements.N && explicit_degree(b_at) == 3
                        if !any(p_at -> p_at.element == Elements.O && explicit_degree(p_at) == 1)
                            push!(beta_atoms, b_at)
                        end
                    end
                end

                if isempty(beta_atoms)
                    nitrogen_count = count(
                        p_at -> p_at.element == Elements.N && explicit_degree(p_at) == 3 && 
                            ( (explicit_valence(p_at) == 4 && is_aromatic(p_at)) ||
                            (explicit_valence(p_at) == 3 && !is_aromatic(p_at)) )
                            , neighbors(at)
                    )
                    # for p_at in (partner(bond, at) for bond in bonds(at))
                    #     if p_at.element == Elements.N && explicit_degree(p_at) == 3
                    #         if ( (explicit_valence(p_at) == 4 && is_aromatic(p_at)) ||
                    #              (explicit_valence(p_at) == 3 && !is_aromatic(p_at)) )
                    #             nitrogen_count += 1
                    #         end

                    #     end
                    # end
                    nitrogen_count >= 2 && (at.atom_type = "80";return)    # Aromatic carbon between N's in imidazolium (CIM+)

                end

                if isempty(beta_atoms) && isempty(alpha_atoms)
                    if at.element == Elements.C
                        if any(p_at -> p_at.element == Elements.C && is_aromatic(p_at) && is_in_ring_6(p_at)
                            , neighbors(at)
                        )
                            at.atom_type = "37";return  #correct atom type for c in c60 (all atoms symmetric)
                                # there is no S:, O:, or N:
                                # this is the case for anions with only carbon and nitrogen in the ring
                                
                        else
                            at.atom_type = "78";return  # General carbon in 5-membered aromatic ring (C5)
                        end
                    elseif at.element == Elements.N
                        at.atom_type = explicit_degree(at) == 3 ? "39" : "76"; return
                        # 39: Aromatic 5 ring nitrogen with pi lone pair (NPYL)
                        # 76: Nitrogen in 5-ring aromatic anion (N5M)
                    end
                end

                if length(alpha_atoms) == 2 &&
                        at.element == Elements.C && in_same_ring(alpha_atoms...) &&
                        alpha_atoms[1].element == Elements.N && alpha_atoms[2].element == Elements.N &&
                        explicit_degree(alpha_atoms[1]) == 3 && explicit_degree(alpha_atoms[2]) == 3;
                    at.atom_type = "80"; return # Aromatic carbon between N's in imidazolium (CIM+)
                end

                if !isempty(alpha_atoms) && isempty(beta_atoms)
                    if      at.element == Elements.C; at.atom_type = "63"; return # Aromatic 5-ring C, alpha to N:, O:, or S: (C5A)
                    elseif  at.element == Elements.N
                        at.atom_type = explicit_degree(at) == 3 ? "81" : "65"; return
                        # 81 : Posivite nitrogen in 5-ring alpha position (N5A+)
                        # 65 : Aromatic 5-ring N, alpha to N:, O:, or S: (N5A)
                    end  
                end

                if isempty(alpha_atoms) && !isempty(beta_atoms)
                    if      at.element == Elements.C; a.atom_type = "64"; return # Aromatic 5-ring C, beta to N:, O:, or S: (C5B)
                    elseif  at.element == Elements.N
                        at.atom_type = explicit_degree(at) == 3 ? "81" : "66"; return
                        # 81 : Posivite nitrogen in 5-ring alpha position (N5B+)
                        # 66 : Aromatic 5-ring N, alpha to N:, O:, or S: (NSB)
                    end  
                end

                if !isempty(alpha_atoms) && !isempty(beta_atoms)
                    if any(pair-> !in_same_ring(pair...), Iterators.product(alpha_atoms, beta_atoms))
                        if      at.element == Elements.C; at.atom_type = "78"; return     # General carbon   in 5-membered aromatic ring (C5)
                        elseif  at.element == Elements.N; at.atom_type = "79"; return # General nitrogen in 5-membered aromatic ring (N5)
                        end 
                    end
                end

                if any(a_at-> a_at.element in (Elements.S, Elements.O), alpha_atoms)
                    if      at.element == Elements.C; at.atom_type = "63";return     # Aromatic 5-ring C, alpha to N:, O:, or S: (C5A)
                    elseif  at.element == Elements.N; at.atom_type = "65";return # Aromatic 5-ring N, alpha to N:, O:, or S: (N5A)
                    end
                end

                if any(b_at-> b_at.element in (Elements.S, Elements.O), beta_atoms)
                    if      at.element == Elements.C; at.atom_type = "64";return     # Aromatic 5-ring C, alpha to N:, O:, or S: (C5A)
                    elseif  at.element == Elements.N; at.atom_type = "66";return     # Aromatic 5-ring N, alpha to N:, O:, or S: (N5A)
                    end
                end

                if      at.element == Elements.C; at.atom_type = "78";return     # General carbon   in 5-membered aromatic ring (C5)
                elseif  at.element == Elements.N; at.atom_type = "79";return     # General nitrogen in 5-membered aromatic ring (N5)
                end
            end

        end

        if is_in_ring_6(at)
            if      at.element == Elements.C;    at.atom_type = "37";return
            elseif  at.element == Elements.N
                if any(p_at -> p_at.element == Elements.O && explicit_degree(p_at) == 1, neighbors(at))
                    at.atom_type = "69"; return # Pyridinium N-oxide nitrogen (NPOX)
                end
                at.atom_type = explicit_degree(at) == 3 ? "58" : "38"; return
                # 58: Aromatic nitrogen in pyridinium (NPD+)
                # 38: Aromatic nitrogen with sigma lone pair (NPYD)     
            end
        end

    end

    if at.element == Elements.H
        for p_at in neighbors(at)
            if p_at.element in (Elements.C, Elements.Si)
                at.atom_type = "5"; return # Hydrogen attatched to carbon (HC) or silicon (HSI)
            end

            if at.element == Elements.O
                if explicit_valence(p_at) == 3
                    at.atom_type = explicit_degree(p_at) == 3 ?  "50" : "52"; return
                    # 50: Hydrogen on oxonium oxygen (HO+)
                    # 52: Hydrogen on oxenium oxygen (HO=+)
                end
    
                for pp_at in neighbors(at)
                    if pp_at.element == Elements.C
                        if is_aromatic(pp_at); at.atom_type = "29";return end    # phenol
                        for ppp_at in neighbors(at)
                            ppp_at.idx == p_at.idx && continue
                            if !(is_aromatic(ppp_at) && is_aromatic(pp_at)) &&
                                    bonded(ppp_at, pp_at).order == BondOrder.Double
                                if      ppp_at.element == Elements.O
                                    at.atom_type = "24";return   # Hydroxyl hydrogen in carboxylic acids (HOCO)
                                elseif  ppp_at.element in (Elements.C, Elements.N)
                                    at.atom_type = "29";return   # Enolic or phenolic hydroxyl hydrogen,
                                                # Hydroxyl hydrogen in HO-C=N moiety (HOCC, HOCN)
                                end
                            end
                        end
                    end
    
                    if pp_at.element == Elements.P; at.atom_type = "24";return end # Hydroxyl hydrogen in H-O-P moiety (HOP)
                    if pp_at.element == Elements.S; at.atom_type = "33";return end # Hydrogen on oxygen attached to sulfur (HOS)
    
                end
    
                hydrogen_count = count(pp_at -> pp_at.element == Elements.H, neighbors(p_at))
                if hydrogen_count == 2; at.atom_type = "32";return end #  Hydroxyl hydrogen in water (HOH)
                at.atom_type = "21";return # Hydroxyl hydrogen in alcohols, Generic hydroxyl hydrogen (HOR, HO)
    
            end

            if p_at.element == Elements.N
                p_at.atom_type == "0" && assign_atom_type_decision_tree(p_at, all_rings, rings_5, rings_6, rings_aro)
                if      p_at.atom_type == "81"; at.atom_type = "36"; return # Hydrogen on imidazolium nitrogen (HIM+)
                elseif  p_at.atom_type == "68"; at.atom_type = "23"; return # Hydrogen on N in N-oxide (HNOX)
                elseif  p_at.atom_type == "67"; at.atom_type = "23"; return # Hydrogen on N in N-oxide (HNOX)
                elseif  p_at.atom_type == "62"; at.atom_type = "23"; return # Generic hydrogen on sp3 nitrogen, e.g., in amines (HNR)
                elseif  p_at.atom_type == "56"; at.atom_type = "36"; return # Hydrogen on guanimdinium nitrogen (HGD+)
                elseif  p_at.atom_type == "65"; at.atom_type = "36"; return # Hydrogen on amidinium nitrogen (HNN+)
                elseif  p_at.atom_type == "43"; at.atom_type = "28"; return # Hydrogen on NSO, NSO2, or NSO3 nitrogen, Hydrogen on N triply bonded to C (HNSO, HNC%)
                elseif  p_at.atom_type == "39"; at.atom_type = "23"; return # Hydrogen on nitrogen in pyrrole (HPYL)
                elseif  p_at.atom_type ==  "8"; at.atom_type = "23"; return # Generic hydrogen on sp3 nitrogen, e.g., in amines, Hydrogen on nitrogen in ammonia (HNR, H3N)
                end

                if explicit_valence(p_at) == 4
                    explicit_degree(p_at) == 2 ? (at.atom_type = "28") : (at.atom_type = "36")
                    return
                    # 28: hydrogen on N triply bonded to C (HNC%)
                    # 36: hydrogen on pyridinium nitrogen, Hydrogen on protonated imine nitrogen (HPD+, HNC+)
                end

                if explicit_degree(p_at) == 2
                    for pp_at in filter(pp_at -> pp_at.element != Elements.H &&
                                                        !(is_aromatic(pp_at) && is_aromatic(p_at) &&
                                                        bonded(pp_at, p_at).order == BondOrder.Double), neighbors(p_at))
                        if pp_at.element in (Elements.C, Elements.N)
                            at.atom_type = "27"; return # Hydrogen on imine nitrogen, Hydrogen on azo nitrogen (HN=C, HN=N)
                        else
                            at.atom_type = "28"; return # Generic hydrogen on sp2 nitrogen (HSP2)
                        end
                    end                                                
                end

                #l. 332
                for pp_at in filter(pp_at -> pp_at.element != Elements.H, collect(neighbors(p_at)))
                    if pp_at.element == Elements.C
                        # deloc. lp pair
                        if is_aromatic(pp_at); at.atom_type = "28"; return end

                        # l. 341
                        if any(ppp_at -> ppp_at.idx != pp_at.idx && 
                                    bonded_aromatic(ppp_at, pp_at) && 
                                    bonded(ppp_at, pp_at).order == BondOrder.Double &&
                                    ppp_at.element in (Elements.C, Elements.N, Elements.O, Elements.S)
                                    , neighbors(pp_at))
                                at.atom_type = "28"; return
                                # Hydrogen on amide nitrogen, Hydrogen on thioamide nitrogen,
                                # Hydrogen on enamine nitrogen, Hydrogen in H-N-C=N moiety (HNCO, HNCS, HNCC, HNCN)
                        end
                    end

                    #l.354
                    if pp_at.element == Elements.N
                        if any(ppp_at -> ppp_at.idx != pp_at.idx && 
                                    bonded_aromatic(ppp_at, pp_at) && 
                                    bonded(ppp_at, pp_at).order == BondOrder.Double &&
                                    ppp_at.element in (Elements.C, Elements.N)
                                    , neighbors(pp_at))
                                at.atom_type = "28"; return
                                # Hydrogen in H-N-N=C moiety, Hydrogen in H-N-N=N moiety (HNNC, HNNN)
                        end
                    end

                    #l.367
                    if pp_at.element == Elements.S
                        if any(ppp_at -> ppp_at.idx != pp_at.idx && 
                                    explicit_degree(ppp_at) == 1 &&
                                    ppp_at.element == Elements.O
                                    , neighbors(pp_at))
                                at.atom_type = "28"; return
                                # Hydrogen on NSO, NSO2 or NSO3 nitrogen (HNSO)
                        end
                    end
                end

                # Generic hydrogen on sp3 nitrogen e.g., in amines,
                # Hydrogen on nitrogen in pyrrole, Hydrogen in ammonia,
                # Hydrogen on N in N-oxide (HNR, HPYL, H3N, HNOX)
                at.atom_type = "23"; return
            end

            # Hydrogen attached to sulfur, Hydrogen attached to >S= sulfur doubly bonded to N,
            # Hydrogen attached to phosphorus (HS, HS=N, HP)
            # if p_at.elemnt in (Elements.S, Elements.P); at.atom_type = "71"; return end
        end     
    end

    # Lithium
    if at.element == Elements.Li
        explicit_degree(at) == 0 && (at.atom_type = "92"; return) # Lithium cation (LI+)
    end

    # Cabron
    if at.element == Elements.C
        if explicit_degree(at) == 4
            if any(mol->any(ring -> length(ring) == 3 && at in ring, mol), all_rings)
                at.atom_type = "22"; return # Aliphatic carbon in 3-membered ring (CR3R)
            end

            if any(mol->any(ring -> length(ring) == 4 && at in ring, mol), all_rings)
                at.atom_type = "20"; return # Aliphatic carbon in 4-membered ring (CR4R)
            end

            at.atom_type = "1"; return # Alkyl carbon (CR)
        end

        if explicit_degree(at) == 3
            bonded_elem = nothing
            oxygen_count = sulphur_count = 0
            N3count = N2count            = 0
            N3fcharge                    = 0
            for p_at in neighbors(at) 
                
                #ToDo: test this; idk if it can be aromatic and double bonded
                if !bonded_aromatic(p_at, at) && bonded(p_at,at).order == BondOrder.Double
                    bonded_elem = p_at.element
                end
        
                if explicit_degree(p_at) == 1
                    p_at.element == Elements.O && (oxygen_count  += 1)
                    p_at.element == Elements.W && (sulphur_count += 1)
                elseif explicit_degree(p_at) == 3
                    p_at.element == Elements.O && (N3count += 1; N3fcharge += p_at.formal_charge)
                elseif explicit_degree(p_at) == 2 && !bonded_aromatic(p_at, at) && 
                    bonded(p_at, at).order == BondOrder.Double && p_at.element == Elements.N
                        N2count += 1
                end
            end

            if N3count >= 2 && !N2count && !oxygen_count && !sulphur_count && 
                (bonded_elem == Elements.N || 
                    !((bonded_elem == Elements.C) && explicit_degree(at) == 3 && N3fcharge == 1))
                at.atom_type = "57"; return # N3==C--N3; Guanidinium carbon, Carbon in +N=C-N: resonance structures (CGD+, CNN+)            
            end

            if oxygen_count == 2 || sulphur_count == 2
                # O1-?-C-?-O1 or S1-?-C-?-S1
                at.atom_type = "41"; return # Carbon in carboxylate anion, Carbon in thiocarboxylate anion (CO2M, CS2M)
            end

            if any(ring-> at in ring, filter(ring -> length(ring) == 3, all_rings)) && bonded_elem == Elements.C
                    at.atom_type = "30"; return # Olefinic carbon in 4-membered ring (CR4E)
            end

            if bonded_elem in (Elements.N, Elements.O, Elements.P, Elements.S)
                # C==N, C==O, C==P, C==S
                at.atom_type = "3"; return
                # Generic carbonyl carbon, Imine-type carbon, Guanidine carbon,
                # Ketone or aldehyde carbonyl carbon, Amide carbonyl carbon,
                # Carboxylic acid or ester carbonyl carbon, Carbamate carbonyl carbon,
                # Carbonic acid or ester carbonyl carbon, Thioester carbonyl (double
                # bonded to O or S), Thioamide carbon (double bonded to S), Carbon
                # in >C=SO2, Sulfinyl carbon in >C=S=O, Thiocarboxylic acid or ester
                # carbon, Carbon doubly bonded to P (C=O, C=N, CGD, C=OR, C=ON, COO,
                # COON, COOO, C=OS, C=S, C=SN, CSO2, CS=O, CSS, C=P)
            end

            at.atom_type = "2"; return # Vinylic Carbon, Generic sp2 carbon (C=C, CSP2)
        end

        explicit_degree(at) == 2 && (at.atom_type = "4" ; return) # Acetylenic carbon, Allenic caron (CSP, =C=)
        explicit_degree(at) == 1 && (at.atom_type = "60"; return) # Isonitrile carbon (C%-)

    end # Cabron end

    # Nitrogen
    if at.element == Elements.N
        if explicit_degree(at) == 4
            any(p_at -> p_at.element == Elements.O && explicit_degree(p_at) == 1) &&
                (at.atom_type = "68"; return) # sp3-hybridized N-oxide nitrogen (N3OX)
            at.atom_type = "34"; return # Quaternary nitrogen (NR+)
        end

        if explicit_degree(at) == 3
            if explicit_valence(at) >= 4
                bonded_elem = nothing
                oxygen_count = count(p_at -> p_at.element == Element.O && explicit_degree(p_at) == 1,neighbors(at))
                #Todo: test this
                if any(p_at -> p_at.element == Elements.N && !bonded_aromatic(at, p_at) &&
                    bonded(at, p_at).order == BondOrder.Double
                    , neighbors(at))
                    bonded_elem = Elements.N
                end

                nitrogen_count = 0
                for p_at in neighbor(at)
                    (bonded_aromatic(at, p_at) || 
                        !(bonded(p_at, at).order == BondOrder.Double)) && continue
                    nitrogen_count += count(
                        pp_at -> pp_at.element == Elements.N && explicit_degree(pp_at) == 3
                        , neighbor(p_at))
                    
                end
                
                if      oxygen_count == 1; at.atom_type = "67"; return # sp2-hybridized N-oxide nitrogen (N2OX)
                elseif  oxygen_count >= 2; at.atom_type = "45"; return # Nitrogen in nitro group, Nitrogen in nitrate group (NO2, NO3)
                end

                if      oxygen_count == 1; at.atom_type = "54"; return # Iminium nitrogen (N+=C)
                elseif  oxygen_count == 2; at.atom_type = "55"; return # Either nitrogen in N+=C-N: (NCN+)
                elseif  oxygen_count == 3; at.atom_type = "56"; return # Guanidinium nitrogen (NGD+)
                end

                bonded_elem == Elements.N && (at.atom_type = "54"; return) # Positivly charged nitrogen doubly bonded to nitrogen (N+=N) 
            end

            if explicit_valence(at) == 3
                is_NNN_or_NNC = false
                is_sulphon_amide = any(p_at -> p_at.element in (Elements.S, Elements.P) &&
                            count(
                                pp_at -> pp_at.element == Elements.O && explicit_degree(pp_at) == 1
                                , neighbors(p_at)) >= 2
                    , neighbors(at))

                is_amide = any(p_at -> p_at.element == Elements.C &&
                        any(pp_at -> !bonded_aromatic(pp_at, p_at) &&
                            bonded(pp_at, p_at).order == BondOrder.Double &&
                            pp_at.element in (Elements.O, Elements.S)
                        , neighbors(p_at))
                    , neighbors(at))
                
                for p_at in neighbors(at)
                    if p_at.element == Elements.C
                        N2count = N3count = 0
                        oxygen_count = sulphur_count = 0

                        # find some properties in nbrnbr
                        double_bonded_elem = triple_bonded_elem = nothing
                        # following functions return the whole atom, not the element of the atoms
                        double_bonded_elem_1 = filter(pp_at -> 
                            !bonded_aromatic(p_at, pp_at) && bonded(pp_at, p_at).order == BondOrder.Double
                            , collect(neighbors(p_at)))
                        double_bonded_elem_2 = filter(pp_at -> 
                            bonded_aromatic(p_at, pp_at) && pp_at.element in (Elements.C, Elements.N)
                            , collect(neighbors(p_at)))
                        triple_bonded_elem = filter(pp_at -> 
                            !bonded_aromatic(pp_at, p_at) && bonded(pp_at, p_at).order == BondOrder.Triple
                            , collect(neighbors(p_at)))
                        
                        if any(pp_at -> pp_at.element == Elements.N && explicit_degree(pp_at) == 3 && 
                                count(ppp_at -> ppp_at.element == Elements.O , neighbors(pp_at)) < 2
                                , neighbors(p_at))
                            N3count += 1
                        end

                        if any(pp_at -> pp_at.element == Elements.N && explicit_degree(pp_at) == 2 &&
                                (bonded(pp_at, p_at).order == BondOrder.Double || bonded_aromatic(pp_at, p_at))
                                , neighbors(p_at))
                            N2count += 1
                        end

                        oxygen_count  += count(pp_at -> is_aromatic(pp_at) && pp_at.element == Elements.O, neighbors(p_at))
                        sulphur_count += count(pp_at -> is_aromatic(pp_at) && pp_at.element == Elements.S, neighbors(p_at))


                        double_bonded_elem_temp = vcat(double_bonded_elem_1, double_bonded_elem_2)
                        double_bonded_elem = !isempty(double_bonded_elem_temp) ? last(double_bonded_elem_temp).element : nothing


                        if N3count == 3; at.atom_type = "56"; return end # Guanidinium nitrogen (NGD+)
                        if !is_amide && !sulphon_amide && !oxygen_count && !sulphur_count && is_aromatic(p_at)
                            at.atom_type = "40"; return 
                        end
                        if N3count == 2 && (double_bonded_elem == Elements.N) && !N2count
                            at.atom_type = "55"; return # Either nitrogen in N+=C-N: (NCN+)
                        end
                    end
                    #638-681
                    if p_at.element == Elements.N
                        nitrogen_count = 0
                        for pp_at in neighbors(p_at)
                            if bonded_aromatic(pp_at, p_at) || !(bonded(pp_at, p_at) == BondOrder.Double)
                                continue
                            end

                            oxygen_count   = count(pp_at -> pp_at.element == Element.O, neighbors(p_at))
                            sulphur_count  = count(pp_at -> pp_at.element == Element.S, neighbors(p_at))
                            nitrogen_count = count(pp_at -> pp_at.element == Element.N, neighbors(p_at))

                            if !oxygen_count && !sulphur_count && nitrogen_count == 1
                                bonded_to_aromatic_c = any(p_at2 -> is_aromatic(p_at2) && p_at2.element == Elements.C && is_in_ring_6(p_at2), neighbors(at))
                                is_NNN_or_NNC = !bonded_to_aromatic_c
                            end
                        end
                    end
                end

                if is_sulphon_amide; at.atom_type = "43"; return end # Sulfonamide nitrogen (NSO2, NSO3)
                if is_amide; at.atom_type = "10"; return end # Amide nitrogen, Thioamide nitrogen (NC=O, NC=S)

                triple_bonded_elem = !isempty(triple_bonded_elem) ? last(triple_bonded_elem).element : nothing
                if double_bonded_elem in (Elements.C, Elements.N, Elements.P) || triple_bonded_elem == Elements.C
                    at.atom_type = "40"; return
                    # Enamine or aniline nitrogen (deloc. lp), Nitrogen in N-C=N with deloc. lp,
                    # Nitrogen in N-C=N with deloc. lp, Nitrogen attached to C-C triple bond
                    # (NC=C, NC=N, NC=P, NC%C)
                end

                if triple_bonded_elem == Elements.N; at.atom_type = "43"; return end # Nitrogen attached to cyano group (NC%N)

                # Nitrogen in N-N=C moiety with deloc. lp
                # Nitrogen in N-N=N moiety with deloc. lp (NN=C, NN=N)
                if is_NNN_or_NNC; at.atom_type = "10"; return end 

                at.atom_type = "8"; return # Amine nitrogen (NR)
                


            end # explicit valence 3
        end # explicit degree 3

        if explicit_degree(at) == 2
            if explicit_valence(at) >= 4
                is_partner_carbon = any(p_at -> p_at.element == Elements.C, neighbors(at))
                is_triple_bond    = any(p_at -> !bonded_aromatic(p_at, at) && bonded(p_at, at).order == BondOrder.Triple, neighbors(at))

                is_partner_carbon && is_triple_bond && (at.atom_type = "61"; return) # Isonitrile nitrogen (NR%)
                at.atom_type = "53"; return # Central nitrogen in C=N=N or N=N=N (=N=)
            end

            if explicit_valence(at) == 3
                if any(p_at -> p_at.element == Elements.O && !bonded_aromatic(p_at, at) &&
                        bonded(p_at, at).order == BondOrder.Double && explicit_degree(p_at) == 1
                        , neighbors(at))
                    at.atom_type = "46"; return #  Nitrogen in nitroso group (N=O)
                end

                if any(p_at -> p_at.element in (Elements.C, Elements.N) && 
                        !bonded_aromatic(p_at, at) && bonded(p_at, at).order == BondOrder.Double
                        , neighbors(at))
                    at.atom_type = "9"; return #  Iminie nitrogen, Azo-group nitrogen (N=C, N=N)
                end

                if any(p_at -> p_at.element == Elements.S &&
                    count(pp_at -> pp_at.element == Elements.O &&
                        explicit_degree(pp_at) == 1
                        , neighbors(p_at)
                    ) >= 2
                )
                    at.atom_type = "43"; return # Sulfonamide nitrogen (NSO2, NSO3)
                end
                    

            end

            if explicit_valence(at) == 2
                if any(p_at -> p_at.element == Elements.S &&
                    count(pp_at -> pp_at.element == Elements.O &&
                            explicit_degree(pp_at) == 1
                            , neighbors(p_at)
                        ) == 1
                , neighbors(at))
                    at.atom_type = "48"; return # Divalent nitrogen replacing monovalent O in SO2 group (NSO)
                end
            end
        end

        if explicit_degree(at) == 1
            for p_at in neighbors(at)
                if !bonded_aromatic(p_at, at) && bonded(p_at, at).order == BondOrder.Triple
                    if any(pp_at -> at.idx != pp_at.idx && 
                            !(p_at.element == Elements.N && pp_at.element == Elements.N 
                                && explicit_degree(pp_at) ==2)
                            ,neighbors(p_at)
                        )
                        at.atom_type = "42"; return # Triply bonded nitrogen (NSP)
                    end
                end

                if p_at.element == Elements.N && explicit_degree(p_at) == 2
                    at.atom_type = "47"; return # Terminal nitrogen in azido group (NAZT)
                end
            end
        end

        at.atom_type = "8"; return # generic amine nitrogen
    end # Nitrogen

    if at.element == Elements.O
        if explicit_degree(at) == 3; at.atomy_type = "49"; return end # Oxonium oxygen (O+)
        if explicit_degree(at) == 2
            hydrogen_count = count(p_at -> p_at.element == Elements.H, neighbors(at))
            if hydrogen_count == 2; at.atom_type = "70"; return end # Oxygen in water (OH2)
            explicit_valence(at) == 3 && (at.atomy_type = "51"; return) # Oxenium oxygen (O=+)
            at.atom_type = "6"; return 
            # Generic divalent oxygen, Ether oxygen, Carboxylic acid or ester oxygen,
            # Enolic or phenolic oxygen, Oxygen in -O-C=N- moiety, Divalent oxygen in
            # thioacid or ester, Divalent nitrate "ether" oxygen, Divalent oxygen in
            # sulfate group, Divalent oxygen in sulfite group, One of two divalent
            # oxygens attached to sulfur, Divalent oxygen in R(RO)S=O, Other divalent
            # oxygen attached to sulfur, Divalent oxygen in phosphate group, Divalent
            # oxygen in phosphite group, Divalent oxygen (one of two oxygens attached
            # to P), Other divalent oxygen (-O-, OR, OC=O, OC=C, OC=N, OC=S, ONO2,
            # ON=O, OSO3, OSO2, OSO, OS=O, -OS, OPO3, OPO2, OPO, -OP)
        end
        if explicit_degree(at) == 1
            for p_at in neighbors(at)
                oxygen_count  = 0
                sulphur_count = 0
                if p_at.element in (Elements.C, Elements.N)
                    oxygen_count  = count(pp_at -> pp_at.element == Elements.O && explicit_degree(pp_at) == 1, neighbors(p_at))
                    sulphur_count = count(pp_at -> pp_at.element == Elements.S && explicit_degree(pp_at) == 1, neighbors(p_at))
                end

                if p_at.element == Elements.H; at.atom_type = "35"; return end # O---H
                if p_at.element == Elements.C
                    if oxygen_count == 2; at.atom_type = "32"; return end # Oxygen in carboxylate group (O2CM)
                    # O--C
                    # Oxide oxygen on sp3 carbon, Oxide oxygen on sp2 carbon (OM, OM2)
                    if !bonded_aromatic(p_at, at) && bonded(p_at, at).order == BondOrder.Single; at.atom_type = "35"; return
                    # O==C
                    # Generic carbonyl oxygen, Carbonyl oxygen in amides,
                    # Carbonyl oxygen in aldehydes and ketones, Carbonyl
                    # oxygen in acids or esters (O=C, O=CN, O=CR, O=CO)
                    else at.atom_type = "7"; return
                    end
                end

                if p_at.element == Element.N
                    # O-?-N-?-O
                    # Oxygen in nitro group, Nitro-group oxygen in nitrate,
                    # Nitrate anion oxygen (O2N, O2NO, O3N)
                    if oxygen_count >= 2; at.atom_type = "32"; return end
                    if !bonded_aromatic(p_at, at) && bonded(p_at, at).order == BondOrder.Single
                        if explicit_degree(p_at) == 2 || explicit_valence(p_at) == 3
                            at.atom_type = "35"; return # O(-)--N
                        else
                            at.atom_type = "32"; return # Oxygen in N-oxides (ONX)
                        end
                    else
                        #is this supposed to be aromatic or something else?
                        ndab = count(b-> b.order in (BondOrder.Double, BondOrder.Aromatic), bonds(p_at))
                        if ndab + explicit_degree(p_at) == 5; at.atom_type = "32"; return # O--N; Oxygen in N-oxides (ONX)
                        else at.atom_type = "7"; return # O==N; Nitroso oxygen (O=N)
                        end
                    end

                end

                if p_at.element == Element.S
                    if sulphur_count == 1; at.atom_type = "32"; return end # O--S
                    # Single terminal oxygen on sulfur, One of 2 terminal O's on sulfur,
                    # One of 3 terminal O's on sulfur, Terminal O in sulfate anion,
                    # (O-S, O2S, O3S, O4S)
                    if !bonded_aromatic(p_at, at) && bonded(p_at,at).order == BondOrder.Single
                        # O--S
                        at.atom_type =  "32"; return # Single terminal oxygen on sulfur, One of 2 terminal O's on sulfur,
                        # One of 3 terminal O's on sulfur, Terminal O in sulfate anion,
                        # (O-S, O2S, O3S, O4S)
                    else
                        # O==S
                        is_sulfoxide = true
                        oxy_sulphur_bonds = count(pp_at -> pp_at.idx != at.idx && pp_at.element == Elements.O, neighbors(p_at))
                        if any(
                            at.idx != pp_at.idx &&
                            (!bonded_aromatic(pp_at, p_at) && bonded(pp_at, p_at).order == BondOrder.Double && pp_at.element == Elements.C && oxy_sulphur_bonds == 1) ||
                            ((pp_at.element == Elements.O && explicit_degree(pp_at) == 1) || (pp_at.element == Elements.N && explicit_degree(pp_at) == 2))
                            , neighbors(p_at)
                        )
                            is_sulfoxide = false
                        end
                        
                        if is_sulfoxide; at.atom_type = "7" ; return # Doubly bonded sulfoxide oxygen (O=S)
                        else             at.atom_type = "32"; return # (O2S, O2S=C, O3S, O4S)
                        end

                        at.atom_type = "32": return # Oxygen in phosphine oxide, One of 2 terminal O's on sulfur,
                        # One of 3 terminal O's on sulfur, One of 4 terminal O's on sulfur,
                        # Oxygen in perchlorate anion (OP, O2P, O3P, O4P, O4Cl)
                    end
                end
            end
        end
    end # Oxygen

    if at.element == Elements.F
        if explicit_degree(at) == 1; at.atom_type = "11"; return end # Fluorine (F)
        if explicit_degree(at) == 0; at.atom_type = "89"; return end # Fluoride anion (F-)
    end

    if at.element == Elements.Na
        at.atom_type = "93"; return # Sodium cation (NA+)
    end

    if at.element == Elements.Mg
        at.atom_type = "99"; return # Dipositive magnesium cation (MG+2)
    end

    if at.element == Elements.Si
        at.atom_type = "19"; return # Silicon (SI)
    end 

    if at.element == Elements.P
        if explicit_degree(at) == 4; at.atom_type = "25"; return end # Phosphate group phosphorus, Phosphorus with 3 attached oxygens,
        # Phosphorus with 2 attached oxygens, Phosphine oxide phosphorus,
        # General tetracoordinate phosphorus (PO4, PO3, PO2, PO, PTET)
        if explicit_degree(at) == 3; at.atom_type = "26"; return end # Phosphorus in phosphines (P)
        if explicit_degree(at) == 2; at.atom_type = "75"; return end # Phosphorus doubly bonded to C (-P=C)
    end

    if at.element == Elements.S
        if explicit_degree(at) == 4; at.atom_type = "18"; return end # Sulfone sulfur, Sulfonamide sulfur, Sulfonate group sulfur,
        # Sulfate group sulfur, Sulfur in nitrogen analog of sulfone
        # (SO2, SO2N, SO3, SO4, SNO)
        if explicit_degree(at) == 3
            double_bond_to_carbon = any(p_at -> !bonded_aromatic(p_at, at) && bonded(p_at, at).order == BondOrder.Double && p_at.element == Elements.C, neighbors(at))
            oxygen_count  = count(p_at -> !bonded_aromatic(p_at, at) && bonded(p_at, at).order == BondOrder.Double && explicit_degree(p_at) == 1 && p_at.element == Elements.O, neighbors(at))
            sulphur_count = count(p_at -> !bonded_aromatic(p_at, at) && bonded(p_at, at).order == BondOrder.Double && explicit_degree(p_at) == 1 && p_at.element == Elements.S, neighbors(at))

            if oxygen_count == 2
                if double_bond_to_carbon
                    at.atom_type = "18"; return # Sulfone sulfur, doubly bonded to carbon (=SO2)
                end
                at.atom_type = "73"; return # Sulfur in anionic sulfinate group (SO2M)
            end

            if oxygen_count && sulphur_count; at.atom_type = "73"; return end # Tricoordinate sulfur in anionic thiosulfinate group (SSOM)
            at.atom_type = "17"; return # Sulfur doubly bonded to carbon, Sulfoxide sulfur (S=C, S=O)

        end

        if explicit_degree(at) == 2
           double_bonded_to_oxygen = any(p_at -> p_at.element == Elements.O && !bonded_aromatic(p_at, at) && bonded(p_at, at).order == BondOrder.Double, neighbors(at)) 
           if double_bonded_to_oxygen; at.atom_type = "74"; return end# Sulfinyl sulfur, e.g., in C=S=O (=S=O)
           at.atom_type = "15"; return # Thiol, sulfide, or disulfide sulfor (S)
        end

        if explicit_degree(at) == 1
            sulphur_count  = length(filter((p_at, pp_at) -> pp_at.element == Elements.S, Iterators.product(neighbors(at), (neighbors(p_at) for p_at in neighbors_(at)), 2)))
            double_bond_to = filter(p_at -> !bonded_aromatic(p_at, at) && bonded(p_at, at).order == BondOrder.Double, neighbors(at)) 
            double_bond_to = !isempty(double_bond_to) ? last(double_bond_to).element : nothing

            if double_bond_to == Elements.C && sulphur_count != 2
                at.atom_type = "16"; return # Sulfur doubly bonded to carbon (S=C)
            end 

            at.atom_type = "72"; return # Terminal sulfur bonded to P, Anionic terminal sulfur,
                                        # Terminal sulfur in thiosulfinate group (S-P, SM, SSMO)
        end
    end

    if at.element == Elements.Cl
        if explicit_degree(at) == 4
            oxygen_count = count(p_at -> p_at.element == Elements.O, neighbors(at))
            if oxygen_count == 4; at.atom_type = "77"; return end # Perchlorate anion chlorine (CLO4)
        end

        if explicit_degree(at) == 1
            at.atom_type = "12"; return # Chlorine (CL)
        end

        if explicit_degree(at) == 0
            at.atom_type = "90"; return # Chloride anion (CL-)
        end
    end

    if at.element == Elements.K
        at.atom_type = "94"; return # Potasium cation (K+)
    end

    if at.element == Elements.Ca
        if explicit_degree(at) == 0
            at.atom_type = "96"; return # Dipositive calcium cation (CA+2)
        end
    end

    if at.element == Elements.Fe
        if at.formal_charge == 2
            at.atom_type = "87"; return # Dipositive iron (FE+2)
        else
            at.atom_type = "88"; return # Tripositive iron (FE+3)
        end
    end

    if at.element == Elements.Cu
        if at.formal_charge == 1
            at.atom_type = "97"; return # Monopositive copper cation (CU+1)
        else
            at.atom_type = "98"; return # Dipositive copper cation (CU+2)
        end
    end

    if at.element == Elements.Zn
        at.atom_type = "95"; return #Dipositive zinc cation (ZN+2)
    end

    if at.element == Elements.Br
        if explicit_degree(at) == 1
            at.atom_type = "13"; return # Bromine (BR)
        end

        if explicit_degree(at) == 0
            at.atom_type = "91"; return # Bromide anion (BR-)
        end
    end

    if at.element == Elements.I
        if explicit_degree(at) == 1
            at.atom_type = "14"; return # Iodine (I)
        end
    end
end


function assign_to(mol, elem_to_params)
    #set some default type for all non-hydro atoms
    # untyped atoms will have atom_type = "0" and name = "Any"
    atoms_ = atoms_df(mol)
    atoms_.atom_type .= "0"         #default value for atom typer
    atoms_.name .= "Any"            #default name
    unassigned_atoms = Dict{String3, Set{Int}}()
    map(eachrow(atoms_)) do at
        push!(get!(unassigned_atoms, string(at.element), Set{Int}()), at.idx)
    end

    # each pair in elem_to_params amps an element to all SMARTS rules for that elem
    # assign the type of an atom if it is caught by one of the rules and
    # if it doesnt already have a type assigned
    for (elem, u_ats) in unassigned_atoms
        for tup in reverse(elem_to_params[elem])
            length(u_ats) == 0 && break
            m = @pipe SMARTSQuery(String(tup.rule)) |> match(_, mol) |> 
                atoms.(_) |>  Iterators.flatten(_)
            isempty(m) && continue
            untyped_atoms = (at for at in m if at.name == "Any" && string(at.element) == elem)
            
            # ignoring assignSpecificValues in atomTyper.C line 216
            map(untyped_atoms) do at
                at.name = tup.symbol
                at.atom_type = tup.type
                delete!(u_ats, at.idx)
            end
        end
    end

end


# `b` corresponds to the "BondStretch" file's column `b`, analogously for c
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

struct BondData{T<:Real}
    kb_normal::T
    r0_normal::T
    standard_bond_exists::Bool
                
    kb_sbmb::T
    r0_sbmb::T
    sbmb_exists::Bool
    empirical::Bool
end

#struct BondData?



function assign_bond_types_mmff94(ac::AbstractAtomContainer{T}, 
        mmff94_params,
        aromatic_rings) where T <: Real
    
    non_hydro_bonds = non_hydrogen_bonds(ac)

    ###### make stretch params
    stretch_params = Dict{Stretch_Idx, BondData{T}}()

    map(eachrow(mmff94_params.sections["BondStretch"].data)) do row
        if Bool(row.a)
            stretch_params[SI(row.b, row.c)] = BondData(T(0), T(0), false, T(row.kb), T(row.r0), true, false)
        else
            stretch_params[SI(row.b, row.c)] = BondData(T(row.kb), T(row.r0), true, T(0), T(0), false, false)
        end
    end

    #######
    map(
        bond -> assign_bond_type(
                        bond,
                        stretch_params,
                        aromatic_rings,
                        mmff94_params.sections["AtomTypeProperties"].data
                    ),
        non_hydro_bonds
    )

    nothing
end


# assignMMFF94BondType, MMFF94.C l428
function assign_bond_type(bond, stretch_params, aromatic_rings, atom_types_df)



    # map a1.type + a2.type -> the entry in stretch_params
    a1_at = atom_by_idx(parent(bond), bond.a1)
    a2_at = atom_by_idx(parent(bond), bond.a2)
    bond_data = get(() -> missing, stretch_params, SI(a1_at.atom_type, a2_at.atom_type))
    ismissing(bond_data) && return false

    is_sbmb::Bool = bond_data.sbmb_exists && !bond_data.standard_bond_exists
    if bond_data.sbmb_exists == bond_data.standard_bond_exists
        if bond.order == BondOrder.Single

            #test this
            #potential problem: atom as from atoms() vs as from atoms_df
            at1_aromatic = a1_at in Iterators.flatten(aromatic_rings)
            at2_aromatic = a2_at in Iterators.flatten(aromatic_rings)

            if !(at1_aromatic || at2_aromatic) &&
                    Bool(atom_types_df[at_idx(a1_at.atom_type), :sbmb]) &&
                    Bool(atom_types_df[at_idx(a2_at.atom_type), :sbmb])
                is_sbmb = true
            elseif at1_aromatic && 
                    at2_aromatic &&
                    !any(ring -> a1_at in ring && a2_at in ring, aromatic_rings) 
                is_sbmb = true
            end
        end
    end
    is_sbmb && set_property!(bond, :MMFF94SBMB, true)
    !is_sbmb && haskey(bond.properties, :MMFF94SBMB) && (delete!(bond.properties, :MMFF94SBMB))

    return true

end


function get_partial_charge_(a,b,c, partial_charges_df, invalid_val, partial_bond_charges_df)
    b == c && return 0
    (b == 0 || c == 0) && return invalid_val

    df = @rsubset partial_charges_df begin
        :a == a
        :b == b
        :c == c
        @kwarg view = true
    end
    
    !isempty(df) && only(df).charge != invalid_val && (return only(df).charge)
    

    #####heuristic value
    p1 = only(partial_bond_charges_df[partial_bond_charges_df.type .== b, :pbci])
    p2 = only(partial_bond_charges_df[partial_bond_charges_df.type .== c, :pbci])

    r = p2 - p1
    r == invalid_val && (@info "No ES parameters bt:$a at1:$b at2:$c."; )
    return r

end

#MMFF94Processors.C line 716
function assign_charges(ac::AbstractAtomContainer{T},
        mmff94_params, aro_rings) where T<:Real

    unassigned_atoms =  []
    invalid_val = 99.0
    # process "Charges"
        # rule_types : whereever :charge is "*". store :symbol
        # types_to_charges: :symbol => :charge
    charge_df = @rsubset mmff94_params.sections["Charges"].data begin
        :charge != "*"
        @kwarg view = true
    end

    rule_charge_df = @rsubset mmff94_params.sections["Charges"].data begin
        :charge == "*"
        @kwarg view = true
    end

    partial_charges_df = mmff94_params.sections["PartialCharges"].data

    # process "PartialBondCharges"
        # phis_ == :fcadj ?
    partial_bond_charges_df = mmff94_params.sections["PartialBondCharges"].data

    assign_partial_charges(ac, charge_df, rule_charge_df, aro_rings)
    charges = Dict{Atom{T}, T}()

    map(at -> setproperty!(at, :InitialCharge, at.charge), atoms(ac))
    for at in atoms(ac)
        nbonds(at) == 0 && continue
        at.atom_type == "0" && continue
        phi = partial_bond_charges_df[parse(Int8, at.atom_type), :fcadj]
        phi == invalid_val && (push!(unassigned_atoms, at))
        phi == 0 && continue
        c = at.charge * phi
        at.charge = (1 - nbonds(at) * phi) * at.charge

        for bond in bonds(at)
            at2 = get_partner(bond, at)
            haskey(charges, at2) ? (charges[at2] += c) : (charges[at2] = c)
        end
    end

    #add all the charges up
    for (at, c) in zip(keys(charges), values(charges))
        at.charge += c
    end
    

    for at1 in atoms(ac)
        charge = copy(at1.charge)
        pcharge = 0

        for bond in bonds(at1)
            at2 = get_partner(bond, at1)
            at1_t = parse(Int, at1.atom_type)
            at2_t = parse(Int, at2.atom_type)
            bt = haskey(bond.properties, :MMFF94SBMB)

            if parse(Int8,at1.atom_type) < parse(Int8,at2.atom_type) 

                pcharge = (-1) * get_partial_charge_(
                    bt,at1_t,at2_t, partial_charges_df,
                    invalid_val, partial_bond_charges_df
                    )
            else
                pcharge = get_partial_charge_(
                    bt,at2_t,at1_t, partial_charges_df,
                    invalid_val, partial_bond_charges_df
                    )
            end
            abs(pcharge) == invalid_val && (push!(unassigned_atoms, at1, at2); continue)
            charge += pcharge

        end

        at1.charge = charge
    end 

    return unassigned_atoms
end


function assign_partial_charges(
        ac::AbstractAtomContainer{T},
        charge_df,
        rule_charge_df,
        aro_rings) where T<:Real


    atoms_df(ac).charge .= 0
    # get all atoms for which the charge can just be read from charge_df
    atoms_simple_charge = @rsubset atoms_df(ac) begin
        !(in(:name, rule_charge_df.symbol))
        @kwarg view = true
    end

    for at in eachrow(atoms_simple_charge)
        idx = findfirst(x -> x == at.name, charge_df.symbol)
        isnothing(idx) && continue
        at.charge = parse(T,charge_df[idx, :charge])
    end

    #these can push atoms into N5Ms and NIMs
    #do N5M find all arom rings which have atom
    
    aro_rings_ = [isempty(ring) ? [] : first(ring) for ring in aro_rings]
    N5Ms = [[ at for at in ring if at.name == "N5M" ] for ring in aro_rings_]
    NIMs = [[ at for at in ring if at.name == "Nim+"] for ring in aro_rings_]

    #do O2S
    atoms_O2S = @rsubset atoms_df(ac) begin
        :name == "O2S"
        @kwarg view = true
    end

    for atom in eachrow(atoms_O2S)
        for bond in bonds(atom)
            partner = get_partner(bond, atom)
            if partner.element == Elements.S && nbonds(partner) == 3
                nr = 0
                map(bonds(partner)) do p_bond
                    !in(p_bond.order, (BondOrder.unknown, BondOrder.Aromatic)) && 
                        (nr += p_bond.order)
                    p_bond.order == BondOrder.Aromatic && 
                        (nr += 5)
                end
                nr == 3 && (atom.charge = -0.5; break)
            end
        end
    end

    #do NR%
    atoms_NR = @rsubset atoms_df(ac) begin
        :name == "NR%"
        @kwarg view = true
    end

    for atom in eachrow(atoms_NR)
        for bond in bonds(atom)
            bond.order == BondOrder.Triple                    && 
                get_partner(bond, atom).element == Elements.N &&
                (atom.charge = 1;)
        end
    end

    #do SM
    atoms_NR = @rsubset atoms_df(ac) begin
        :name == "SM"
        @kwarg view = true
    end

    for atom in eachrow(atoms_NR)
        partner = get_partner(only(bonds(atom)), atom)
        for bond in bonds(partner)
            partner2 = get_partner(bond, partner)
            if partner2.element == Elements.O && nbonds(partner2) == 1 &&
                only(bonds(partner2)).order == BondOrder.Single
                atom.charge = -0.5
            end
        end

        arom.charge = -1
    end

    #handle nim, n5m
    for ring in (r for r in N5Ms if !isempty(r))
        charge_ = -1/length(ring)
        map(at -> at.charge = charge_, ring)
    end

    for ring in (r for r in NIMs if !isempty(r))
        nr_partner = 0
        for at in ring
            for bond in bonds(ring)
                get_partner(bond, at).name in ("NGD+", "NCN+") && (nr_partner +=1) 
            end
        end

        charge_ = 1.0 / (length(ring) + nr_partner)
        map(at -> at.charge = charge_, ring)
    end

end