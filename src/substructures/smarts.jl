using MolecularGraph: 
    QueryMol, 
    smartstomol, substructmatches
export SMARTSQuery

struct SMARTSQuery
    query::String
    query_graph::QueryMol

    function SMARTSQuery(query::String)
        new(query, smartstomol(query))
    end
end

function _to_substructure(name, mol::AbstractMolecule{T}, m; adjacent_bonds=false) where {T<:Real}
    matched_atoms = keys(m)

    filter_atoms(
            n -> getfield(n, :rownumber) ∈ matched_atoms, mol; 
            name=name, adjacent_bonds=adjacent_bonds
    )
end

function Base.match(query::SMARTSQuery, mol::AbstractMolecule{T}; adjacent_bonds=false) where {T<:Real}
    mg_mol = convert(GraphMol{SDFileAtom, SDFileBond}, mol)
    matches = substructmatches(mg_mol, query.query_graph)

    [_to_substructure("$(query.query) on $(mol.name)", mol, m; adjacent_bonds) for m in matches]
end