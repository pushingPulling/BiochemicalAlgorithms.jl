import MolecularGraph: 
    MolGraph, 
    SDFAtom, 
    SDFBond,
    sssr,
    isringatom

export find_sssr, is_ring_atom

function _filter_bonds(ac::AbstractAtomContainer{T}) where {T<:Real}
     # filter out cystein bridges and h bridges (TODO!)
     new_atoms = filter(x -> true, _atoms(ac), view=true)
     new_bonds = filter(b -> !get(b.properties, "DISULPHIDE_BOND", false), 
                        _bonds(ac), view=true)

    convert(MolGraph{SDFAtom, SDFBond}, Substructure(ac.name, ac, new_atoms, new_bonds, ac.properties))
end


function find_sssr(ac::AbstractAtomContainer{T}) where {T<:Real}
    mg = _filter_bonds(ac)

    mg_sssr = sssr(mg)

    map(r->map(a->atom_by_number(ac isa System ? ac : parent_system(ac), collect(keys(mg.vprops))[a]), r), mg_sssr)
end

function is_ring_atom(ac::AbstractAtomContainer{T}) where {T<:Real}
    mg = _filter_bonds(ac)

    isringatom(mg)
end