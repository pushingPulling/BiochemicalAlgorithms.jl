export load_sdfile, write_sdfile

using MolecularGraph: sdfilereader, sdfilewriter

function load_sdfile(fname::String, T=Float32)
    mg_mols = [m for m in sdfilereader(fname)]

    sys = System{T}(fname)

    for mg_mol in mg_mols
        convert(Molecule{T}, mg_mol; system=sys)
    end

    sys
end

function write_sdfile(fname::String, mol::AbstractMolecule)
    mg_mol = convert(MolGraph{SDFAtom, SDFBond}, mol)

    sdfilewriter(fname, [mg_mol])
end

function write_sdfile(fname::String, sys::System)
    mg_mols = [convert(MolGraph{SDFAtom, SDFBond}, m) for m in molecules(sys)]

    sdfilewriter(fname, mg_mols)
end