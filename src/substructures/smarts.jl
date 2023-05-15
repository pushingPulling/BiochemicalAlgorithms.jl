using MolecularGraph: 
    smartstomol, substruct_matches
export SMARTSQuery

struct SMARTSQuery
    query::String
    query_graph::MolGraph

    function SMARTSQuery(query::String)
        new(query, smartstomol(query))
    end
end

function _to_substructure(name, mol, m; adjacent_bonds=false)
    matched_atoms = keys(m)

    filter_atoms(
            :number => n -> n ∈ matched_atoms, mol; 
            name=name, adjacent_bonds=adjacent_bonds
    )
end

function Base.match(query::SMARTSQuery, mol::AbstractMolecule; adjacent_bonds=false)
    mg_mol = convert(MolGraph{SDFAtom, SDFBond}, mol)
    matches = substruct_matches(mg_mol, query.query_graph)

    [_to_substructure("$(query.query) on $(mol.name)", mol, m; adjacent_bonds) for m in matches]
end