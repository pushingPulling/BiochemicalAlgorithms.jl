export kekulizer!

mutable struct AtomInfo
    atom::Atom
    double_bond::Union{Bond, Missing}
    abonds::Vector{Bond}
    partner_idxs::Vector{Int}
    cur_double::Int
    min_double::MaybeInt
    max_double::MaybeInt
    uncharged_double::MaybeInt
    max_valence::MaybeInt

    #original c++ code initializes
    function AtomInfo(atom_::Atom)
        return new(atom_, missing, Bond[], Int[], 0, missing, missing, missing, missing)
    end
end


#new implementation of C++ BALL "Kekulizer::setup"
#the input is an iterable of aromatized rings that are present in `mol`!
#one way to aquire aromatized rings is qsar::aromaticity::aromatize_simple
function kekulizer!(mol, aro_rings)
    # i chose to save the information about the maximul valence for each atom twice:
    # once in the atom itself as a property (in atom.properties) and once in each atom_info

    # a problem that came up: the atom_ifnos that make up 'solutions' need to be deepcopied, 
    # which also mean that in 'apply_solution' we must access the original bonds_df of the system in an ugly way


    #get SSSR & aromatic "systems" - system here is generalization of ring.
    fix_substructures(mol)

    ok::Bool = all(map(fix_aromatic_ring, aro_rings))

    #set unassigned bonds
    unassigned_bonds = filter(bond -> bond.order == BondOrder.Aromatic, bonds(mol))

    #set formal charges
    if ok
        for at in Iterators.flatten(aro_rings)
            nr = 0
            for bond in bonds(at)
                bond.order in (BondOrder.Single, BondOrder.Double, 
                                BondOrder.Triple, BondOrder.Quadruple) &&
                    (nr += Int(bond.order))
                bond.order == BondOrder.Aromatic && (nr += 5)
            end
            at.formal_charge = nr - pop!(at.properties, :kekulizer_max_valence)
        end
    end

    #i should return the unassigned bonds here?
    return unassigned_bonds

end

function fix_substructures(mol)
    #find and process substructures using smarts matcher
        #handle carboxilic acid

    #original query in C++ BALL was [#8;D1]~#6~[#8;D1]
    #but that sequence is just molecular oxygen?
    m = unique(match(SMARTSQuery("[#8;D1]~[#8;D1]"), mol))

    #this query is actually carboxilic acid according to 
    #https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html
    #the smarts result is buggy. the match doesnt capture all oxygens
    m = unique(match(SMARTSQuery("[CX3](=O)[OX2H1]"), mol))


    for match in m
        carbon = filter(at -> at.element == Elements.C, atoms(match))
        oxygens = filter(b -> any(at -> at.element == Element.O, (b.a1, b.a2)), bonds(carbon))

        !any(b -> b.order == BondOrder.Aromatic, bonds(carbon)) && continue

        #one oxygens has only a single bond
        if nbonds(first(oxygens)) == 1
            map(o -> o.order = BondOrder.Single, oxygens)
            first(oxygens).order = BondOrder.Double
        else
            map(o -> o.order = BondOrder.Single, oxygens)
            oxygens[2].order = BondOrder.Double
        end

        #erase atoms? nr_ca++?
    end

    #fix amidine and guanine
    m = unique(match(SMARTSQuery("[#7;D1]~[#6R0]~[#7;D1]"), mol))
    
    for match in m
        carbon = filter(at -> at.element == Elements.C, atoms(match))
        nitrogens = filter(b -> any(at -> at.element == Element.N, (b.a1, b.a2)), bonds(carbon))

        !any(b -> b.order == BondOrder.Aromatic, bonds(carbon)) && continue

        #if bonds are with C and H -> double
        if nbonds(first(nitrogens)) == 2 && 
            count(b -> any(at -> at.element == Elements.H, (b.a1, b.a2)), bonds(first(nitrogens))) == 1

            map(o -> o.order = BondOrder.Single, nitrogens)
            first(nitrogens).order = BondOrder.Double
        else
            map(o -> o.order = BondOrder.Single, nitrogens)
            nitrogens[2].order = BondOrder.Double
        end

        #erase atoms and nr_am_gu up?
    end

    #fix phosphoric acid

    m = unique(match(SMARTSQuery("[P]([#8;D1])([#8;D1])([#8;D1])"), mol))

    for match in m
        carbon = filter(at -> at.element == Elements.C, atoms(match))
        oxygens = filter(b -> any(at -> at.element == Element.O, (b.a1, b.a2)), bonds(carbon))

        !any(b -> b.order == BondOrder.Aromatic, bonds(carbon)) && continue

        #one oxygen has only a single bond
        oxygen_double_bond_idx = oxygens[findfirst(ox -> nbonds(ox) == 1, oxygens)]
        
        map(o -> o.order = BondOrder.Single, oxygens)
        oxygens[oxygen_double_bond_idx].order = BondOrder.Double
        
        #erase oxygens and nr_phos ++?
        
    end

    #in C++ ball the smart matches something like O~O~C
    #but the comment mentions O~N and only sets one bond.
    m = unique(match(SMARTSQuery("[OR0D1\$(O~[CR0])]"), mol))

    for match in m
        oxygen = atoms(match)[findfirst(at -> nbonds(at) == 2, atoms(match))]

        !any(b -> b.order == BondOrder.Aromatic, bonds(oxygen)) && continue

        map(b -> b.order = BondOrder.Double, bonds(oxygen))

        #count O~N connections?
        
    end
end


#original fixAromaticRings in kekulizer.C@274
function fix_aromatic_ring(aro_ring)

    #abort for structures that do not look like a ring
    length(aro_ring) < 3 && return false

    at_infos = Dict{Int, AtomInfo}(at.idx => AtomInfo(at) for at in aro_ring)
    max_valence!(at_infos)


    for at_info in values(at_infos)
        cur_valence = 0
        for bond in non_hydrogen_bonds(at_info.atom)
            if !in(Int(bond.order), 
                (BondOrder.Double, BondOrder.Triple, BondOrder.Quadruple))
                cur_valence += 1
            else
                #double = 2; triple = 3; quadruple =4
                cur_valence += Int(at_info.atom.order)
            end
        end

        at_info.uncharged_double = at_info.max_valence - cur_valence
        if at_info.uncharged_double < 0
            @info "Kekulizer could not 
                calculate max number of needed double bonds for $(at_info.atom)."
            return false
        end
        at_info.max_double = min(1, at_info.uncharged_double + 1)
        at_info.min_double = max(0, at_info.uncharged_double - 1)

        if !in(Int(at_info.atom.element), (6,7))
            at_info.max_double = at_info.uncharged_double
            at_info.min_double = at_info.uncharged_double
        end

        for bond in bonds(at_info.atom)
            bond.order == BondOrder.Aromatic && 
                at_info.atom.idx > get_partner(bond, at_info.atom).idx &&
                (bond.order = BondOrder.Single; push!(at_info.abonds, bond))
        end
        
    end 



    solutions = []
    #using a Ref here so i can use this value inside this func and fix_aromatic_systems
    lowest_penalty = Ref{Int}(typemax(Int))
    cur_penalty = 0
    fix_aromatic_systems(lowest_penalty, cur_penalty, solutions, at_infos, 1)

    if lowest_penalty[] < typemax(Int)
        if lowest_penalty[] == 0
            rows = (atom_info.double_bond for atom_info in only(solutions))
            idxs = [getfield(getfield(row, :_row), :dfrow) for row in rows]
            bonds_df(parent(first(aro_ring)))[idxs, :order] .= BondOrder.Double
        else
            rows = (atom_info.double_bond for atom_info in values(distance_scores(solutions)))
            idxs = [getfield(getfield(row, :_row), :dfrow) for row in rows]
            bonds_df(parent(first(aro_ring)))[idxs, :order] .= BondOrder.Double
        end
    else
        return false
    end

    return true

end



function distance_scores(solutions)
    best_score = typemax(Int)
    best_solution = 0

    length(solutions) == 1 && return only(solutions)

    for solution in solutions

        pos_atoms = []
        neg_atoms = []

        for at_info in values(solution)
            at_info.cur_double < at_info.uncharged_double && (push!(neg_atoms, at_info.atom))
            at_info.cur_double > at_info.uncharged_double && (push!(pos_atoms, at_info.atom))
        end

        score = 0
        for i in eachindex(pos_atoms)
            for j in eachindex(pos_atoms)
                score += 100 / (distance(pos_atoms[i], pos_atoms[j]) + 1)
            end
        end

        for i in eachindex(neg_atoms)
            for j in eachindex(neg_atoms)
                score += 100 / (distance(neg_atoms[i], neg_atoms[j]) + 1)
            end
        end

        for i in eachindex(pos_atoms)
            for j in eachindex(neg_atoms)
                score += (distance(neg_atoms[i], neg_atoms[j])) / 100
            end
        end
        
        score < best_score && (best_score = score; best_solution = solution)
    end

    return values(best_solution)
end


function fix_aromatic_systems(lowest_p, cur_penalty, solutions, at_infos, idx)

    cur_penalty > lowest_p[] && return

    if idx > length(values(at_infos))      #this is faster than collecting all values
        lowest_p[] == 0 && return       #if there is 0 penatly, then no change needed!
        lowest_p[] = min(lowest_p[], cur_penalty)
        cur_penalty < lowest_p[] && (empty!(solutions); lowest_p[] = cur_penalty)
        cur_penalty <= lowest_p[] && (push!(solutions, deepcopy(at_infos)))
        return 
    end


    ai = collect(values(at_infos))[idx]

    if isempty(ai.abonds)
        tap = 0
        ai.cur_double < ai.uncharged_double && (tap = get_penalty(ai, -1))
        cur_penalty += tap
        fix_aromatic_systems(lowest_p, cur_penalty, solutions, at_infos, idx +1)
        cur_penalty -= tap
        return
    end 

    #test if uncharged solutionn is available
    if ai.min_double <= ai.cur_double && ai.uncharged_double == ai.cur_double
        fix_aromatic_systems(lowest_p, cur_penalty, solutions, at_infos, idx +1)
    end

    # try to find solution by adding another double bond
    tap = 0
    if ai.cur_double < ai.max_double
        ai.uncharged_double != 1 && (tap = get_penalty(ai, 1))

        for abond in ai.abonds
            partner = at_infos[get_partner(abond, ai.atom).idx]
            partner.cur_double == partner.max_double && continue

            pap = 0             #partner atom penalty
            partner.cur_double + 1 > partner.uncharged_double && 
                (pap += get_penalty(partner, 1))
            pap += tap
            cur_penalty + pap > lowest_p[] && continue

            ai.cur_double += 1
            partner.cur_double += 1
            ai.double_bond = abond
            partner.double_bond = abond
            cur_penalty += pap

            fix_aromatic_systems(lowest_p, cur_penalty, solutions, at_infos, idx +1)

            ai.cur_double -= 1
            partner.cur_double -= 1
            ai.double_bond = missing
            partner.double_bond = missing
            cur_penalty -= pap

            
        end
    end
    
    if ai.uncharged_double != ai.cur_double
        tap = 0
        ai.cur_double == 0 && (tap = get_penalty(ai, -1))

        cur_penalty += tap
        fix_aromatic_systems(lowest_p, cur_penalty, solutions, at_infos, idx +1)
        cur_penalty -= tap
    end
        
end


function get_penalty(ai, charge; use_formal_charges = true)
    POSITIVE_NITROGEN = 10
    NEGATIVE_NITROGEN = 11
    NEGATIVE_CARBON   = 25
    POSITIVE_CARBON   = 26
    UNEQUAL_CHARGE    = 100

    at = ai.atom
    if use_formal_charges && at.formal_charge != 0
        at.formal_charge != charge && return UNEQUAL_CHARGE
        return 0
    end

    if at.element == Elements.C 
        charge ==  1 && return POSITIVE_CARBON
        charge == -1 && return NEGATIVE_CARBON
    end

    if at.element == Elements.N 
        charge == -1 && return NEGATIVE_NITROGEN
        charge == 1 && any(b -> nbonds(get_partner(b, at)) == 1, bonds(at)) &&
            return POSITIVE_NITROGEN -1
        return POSITIVE_NITROGEN
    end

    return 0

end 


function max_valence!(at_infos)
    max_valence_set = Set{Atom}()

    for at_info in values(at_infos)
        at_info.atom in max_valence_set ? continue : (push!(max_valence_set, at_info.atom))

        if Int(at_info.atom.element) == 6 
            at_info.max_valence = 4
            at_info.atom.properties[:kekulizer_max_valence] = 4
        elseif Int(at_info.atom.element) in (8,16,34,52)
            at_info.max_valence = 2
            at_info.atom.properties[:kekulizer_max_valence] = 2
        elseif Int(at.element) in (7,15,33)
            at_info.max_valence = 3
            at_info.atom.properties[:kekulizer_max_valence] = 3
        end

        if Int(at_info.atom.element) in (7,16)
            for bond in bonds(at_info.at)
                at_info.atom.order == BondOrder.Double && 
                    Int(get_partner(bond, at_info.atom).element) == 8 &&
                    (at_info.atom.max_valence += 2; at_info.atom.properties[:kekulizer_max_valence] += 2)
            end
        end

    end
end 

