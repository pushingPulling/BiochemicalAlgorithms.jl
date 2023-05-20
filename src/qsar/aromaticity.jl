export
    simple_can_be_aromatic_,
    simple_can_be_aromatic_weaker_

# find all rings which contain atom `at`
function rings_(at::Atom{T}, rings) where T
    filter(ring -> (in(at, ring)), rings)
end

function hueckel_succeeds_(ring)
    (count_pi_electrons_(ring) - 2) % 4 == 0
end

function simple_can_be_aromatic_(ring)
    destab = 0
    #for all atoms
    for atom in ring
        single_bonds = 0; double_bonds = 0; aromatic_bonds = 0
        #for each bond adjacent to this atom
        for bond in bonds(atom)
            #if the ring has the partner of the bond
            if get_partner(bond, atom) in ring
                #count s/double/aroma bonds
                bond.order == BondOrder.Single   && (single_bonds   += 1)
                bond.order == BondOrder.Double   && (double_bonds   += 1)
                bond.order == BondOrder.Aromatic && (aromatic_bonds += 1)
            end
        end

        #if element is C
        if atom.element == Elements.C
            if !((double_bonds == 1 && single_bonds > 0) || aromatic_bonds > 1) 
                return false
            end
        else
            atom.element == Elements.S && nbonds(atom) > 2 && return false
            double_bonds == 1 || aromatic_bonds > 1 && (destab += 1)
        end
    end

    return destab < 2
end


function simple_can_be_aromatic_weaker_(ring)
    destab = 0
    #for all atoms
    for atom in ring
        single_bonds = 0; double_bonds = 0; aromatic_bonds = 0
        #for each bond adjacent to this atom
        for bond in bonds(atom)
            bond.order == BondOrder.Single   && (single_bonds   += 1)
            bond.order == BondOrder.Aromatic && (aromatic_bonds += 1)
            if bond.order == BondOrder.Double
                get_property(get_partner(bond, ring), :in_ring) && (double_bonds += 1)
            end


        end

        #if element is C
        if atom.element == Elements.C
            if !((double_bonds == 1 && single_bonds > 0) || aromatic_bonds > 1) 
                return false
            end
        else
            double_bonds == 0 && (destab += 1)
            atom.element == Elements.S && nbonds(atom) > 2 && return false
        end
    end

    return destab < 2
end

function count_pi_electrons_(ring)
    num_pi       = 0; hetero_count   = 0
    single_bonds = 0; double_bonds   = 0
    triple_bonds = 0; aromatic_bonds = 0

    #we know that the ring system has alternating double bonds
    for atom in ring
        # handles charged atoms, tests if the charge is an integer value
        # this is bc BALL knows nothing about formal charges, explicitly
        # i.e. cyclopropyl cation or tropylium cation
        if atom.charge != 0
            atom.charge ==  1 && (num_pi -= 1)
            atom.charge ==  2 && (num_pi -= 2)
            atom.charge ==  3 && (num_pi -= 3)
            atom.charge == -1 && (num_pi += 1)
            atom.charge == -2 && (num_pi += 2)
            atom.charge == -3 && (num_pi += 3)
        end

        #assuming enum "Elements" assigns names to atomic number
        atom.element == 5 && nbonds(atom) > 3 && return 0
        if atom.element in (6, 14, 32, 50)
            single_bonds = 0; double_bonds = 0; triple_bonds = 0
            for bond in bonds(atom)
                bond.order == BondOrder.Double   && (double_count   += 1)
                bond.order == BondOrder.Aromatic && (aromatic_bonds += 1)

                bond.order == BondOrder.Triple && 
                    get_partner(bond, atom).element == Elements.C && 
                    (triple_bonds += 1)
                
            end

            if double_bonds ==1 || triple_bonds == 1 || aromatic_bonds == 2
                num_pi += 1
            else
                return 0
            end

        elseif atom.element in (7,15,33,51) #N, P, As, Sb  - very experimental for As, Sb
            single_bonds = 0; double_bonds = 0; triple_bonds = 0
            
            for bond in bonds(atom)
                bond.order == BondOrder.Single   && (single_bonds   += 1)
                bond.order == BondOrder.Double   && (double_count   += 1)
                bond.order == BondOrder.Aromatic && (aromatic_bonds += 1)
            end

            # if this test is true we have a problem, maybe a P or N
            # which has more than 3 bonds
            double_bonds > 1  || single_bonds > 3 && return 0
            double_bonds == 1 || (aromatic_bonds == 2 && single_bonds == 0) &&
                (num_pi +=1) 
            double_bonds == 0 || (aromatic_bonds == 2 && single_bonds == 1) &&
                (num_pi += 2; hetero_count += 1)

        elseif atom.element in (8,16,34,52)
            num_pi += 2
            hetero_count += 1
        
        # C++ Ball throws an error at this point
        else
            throw("aromaticity.jl/count_pi_electrons_: No pi-electron \
            handle for atom with element $(atom.element)")
        end
    end

    hetero_count > 1 && return 0
    return num_pi

end

#aromatize simple. one ring -> aromatized ring
function aromatize_simple(ring)

    atom_to_rings = Dict{}()
    aromatic_rings = []
    can_be_rings = []
    can_be_aromatic = simple_can_be_aromatic(ring)
    can_be_aromatic_weaker = can_be_aromatic ? false : simple_can_be_aromatic_weaker_(ring)

    
    if can_be_aromatic || can_be_aromatic_weaker
        
        if hueckel_succeeds_(ring)
            #test if the rings contains a sp3 nitrogen
            has_sp2n::Bool = false
            for atom in ring
                if atom.element == Elements.N
                    double_bonds = count(bond -> bond.order == BondOrder.Double, bonds(atom))
                    double_bonds < 1 && (atom_to_rings[atom] = ring; has_sp2n = true)
                end
            end
            ring_container = can_be_aromatic ? aromatic_rings : can_be_rings
            !has_sp2n ? push!(ring_container, ring) : push!(sp2n_rings, ring)   
        end

    end



    #handle the sp2n containing rings - l179
    for atom in keys(atom_to_rings)
        rings = atom_to_rings[atom]
        if length(rings) == 2
            if simple_can_be_aromatic_(rings[1]) && !simple_can_be_aromatic_(rings[2])
                if hueckel_succeeds_(rings[1])
                    push!(aromatic_rings, rings[1])
                else
                    simple_can_be_aromatic_weaker_(rings[2]) &&
                        hueckel_succeeds_(rings[2]) &&
                        push!(can_be_rings, rings[2])
                end
                continue;
            end

            if !simple_can_be_aromatic_(rings[1]) && simple_can_be_aromatic_(rings[2])
                if hueckel_succeeds_(rings[2]) 
                    push!(aromatic_rings, rings[2])
                else    #215
                    simple_can_be_aromatic_(rings[1]) &&    #unlike in codeblock above, this one doesnt check weaker
                        hueckel_succeeds_(rings[1])   &&
                        push!(aromatic_rings, rings[1])
                end
                continue;
            end

            
            if simple_can_be_aromatic_(rings[1]) && simple_can_be_aromatic_(rings[2])
                success = findfirst(r -> hueckel_succeeds_(r), rings)
                push!(aromatic_rings, rings[success])
                continue;
            end

            
            success = findfirst(r -> simple_can_be_aromatic_weaker_ && hueckel_succeeds_(r), rings)
            push!(can_be_rings, rings[success])
            continue;
                
            
        elseif length(rings) == 3
            success = findfirst(r -> simple_can_be_aromatic_(r) && hueckel_succeeds_(r), rings)
            !isnothing(success) && (push!(aromatic_rings, rings[success]); continue;)
            success = findfirst(r -> simple_can_be_aromatic_weaker_(r) && hueckel_succeeds_(r), rings)
            !isnothing(success) && (push!(can_be_rings, rings[success]); continue;)
            
        elseif length(rings) == 1
            simple_can_be_aromatic_(rings[1]) &&
                hueckel_succeeds_(rings[1])   &&
                (push!(aromatic_rings, rings[1]); continue;)
            simple_can_be_aromatic_weaker_(rings[1]) &&
                push!(can_be_rings, rings[1])
        end
    end
    
    # now handle the rings which can be aromatic 
    aromatic_atoms = [at for at in ring for ring in aromatic_rings]
    for ring in can_be_aromatic
        can_be::Bool = true
        for atom in ring
            single_bonds = 0; double_bonds = 0; aromatic_bonds = 0;
            for bond in bonds(atom)
                if get_partner(bond, atom) in ring
                    bond.order == BondOrder.Single   && (single_bonds   += 1)
                    bond.order == BondOrder.Double   && (double_bonds   += 1)
                    bond.order == BondOrder.Aromatic && (aromatic_bonds += 1)
                else
                    get_partner(bond, atom) in aromatic_atoms && (aromatic_bonds += 2)
                end
            end

            atom.element == Elements.C &&
                !((double_bonds == 1 && single_bonds > 0) || aromatic_bonds > 1) &&
                (can_be = false; break;)
        end

        if can_be
            for aro_ring in aromatic_rings      #why iterate over this? its like that in original ball
                if hueckel_succeeds_(ring)
                    push!(aromatic_rings, ring)
                    break
                end
            end
        end
    end


end