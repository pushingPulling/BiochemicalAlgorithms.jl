# read up on
-   Born Oppenheimer
-   CCSD(T)
-   Hartree fock (molecular orbital)

# MOOP plot
at 10k reps, step=10

no moop version:
mean(norm.(x[1] .- x[3])) = 1.9903777722551101
mean(norm.(x[3] .- x[4])) = 1.948130567640203 
mean(norm.(x[4] .- x[1])) = 1.7792434546482523

moop version
mean(norm.(x[1] .- x[3])) = 1.9356008590398346
mean(norm.(x[3] .- x[4])) = 1.9777156568343008
mean(norm.(x[4] .- x[1])) = 1.8868945320245636

# testing
8ce0 + cpx351
atoms = 1432
Stretches: 1269
Bends: 1717
Stretch-Bends: 1716
Torsions: 2002
Out-of-plane Bends: 225
ES interactions: 279215
VDW interactions: 269876

no optims bench - 3 iters, step=10 ,displace=(46.5, 30,74)
10.212 sec

Profile:
ES: 41.2%               -
VDW: 44.3%              - 88% in total
            most of the time in that is read/write

# ML-FF
learn FF params from ab initio data: need only atomic positions(not bonds)


# more references
i n gpu_amber 2 at the start of esction 3 there are more papers about GPU PME for MDS



# MDS 
 Berendsen weak coupling algorithm,43 the Andersen temperature coupling scheme,44 and the Langevin dynamics thermostat



# stuff for the presentation
- insert all the plots
- practice presentation/ change some things
    - justify that FF is lame with MDS being cool
- read parallel computing paper
- read Machine Learning paper



## from papers
# timescales 
neuroscience: signal transmission on microsecs scale, diffusion microsecs
millisecs: folding, 100ms: catalytic turnover, ligand binding

## MMFF94 description

in MMFF94 paper sind die params und die Energy func beschrieben. [paper](https://onlinelibrary.wiley.com/doi/10.1002/(SICI)1096-987X(199604)17:5/6%3C490::AID-JCC1%3E3.0.CO;2-P)

Es sind die atom types gelistet

MMFF94 ist für MDS geeignet, nicht für energy min

## Forcefields in Biochem

`AbstractForceFieldComponent` is for BondAngles, BondBends, torsion etc
`ForceField` enthält ein system, parameter und mehrere `AbstractForceFieldComponent`
in amber.jl wird ein FF erstellt. Das sollte das Ziel sein
Es gibt `compute_forces` und `compute_energy` was man auf die FF einwirkt.

## Cleanups

- rewrite assign_bond_types to not use BondData.
- use grouped dataframes to search stuff?
- return missing instead of -1


# Testing
- when 2 atoms then there are no bends, no stretch-bends; only 1 stretch exists


# nonbonded

- Metagraph?? nonbonded comp 307
- periodic_box
- cell list map

- Br and Zinc get assigned wrong atom types and charges by C++ Smartmatcher
- the matcher matches [BR;v0] when it shouldnt

# FF Worklist

- mark all the units in the structs.
- correctly handle unassigned atoms
    - check each component if they prodice unassigned atoms
    - check mmff94 setup for the same
    - when i dont know atom type, is it linear or not? should i not produce a bend?
- does VdW Switch have any function? i am not using it at all


# Weiteres Vorgehen

- Kekulizer and mmff94 assign charges - how to handle unassigned atoms?
- components implementieren

# Implement diary
-   füge molgraph 0.14.2 hinzu wegen kekulization
    -   interface hat sich geändert, graphs werden jetzt benutzt. hab imprinzip nur am umschreiben (aka einarbeiten) gearbeitet.
-   abgebrochen, zu inkonsistent und das ganze interface hat sich bei Molgraph 0.13 -> 0.14 geändert
- mache ab jetzt mit molgraph@0.13 weiter
- implementiere c++ BALL'S kekulizer, aromatizer
- wtf macht l179 aromaticiyProcessor
    - debug/look in c++ BALL?
- aromatize_simple fertig
- aromatize_simple test: ((Sucks, cause ball.jl  bond detection is incomplete))
- done actually!
- kekulisierer integrieren: done!
- GraphMol SMARTSQuery has a shitty error: it doesnt detect a double bonda as a double bond
- query in question : "[\$([#7]=[#6])]" doesnt work. But "[\$([#7]-[#6])]" works stoopid program (solution: "=" and "-" in aro ring is ":" in smartsparsers)
- remove bonded and geminal interactions

### Full File MM setup, energy and force test
- sdfile_test1.sdf doesnt work on C++- On julia it has partially missing result (See [\$(#7=#6)] issue)
- not possible on all files due to inconsistent amount of torsions on c++'s end

### untested
- assign hetcyc 5 types for anyhing that isnt C or N

### inconsistencies
- in C++ BALL, on descriptors_test.sdf the sum of the amount of torsions per molecule adds up to 330. When using Setup on all molecules at the same time, then it only adds up to 326. Julia BALL reports 330 in both cases
- on descriptors_test C++ BALL reports atom type 2 (C in a double bond) instead of 37 (C in an aromatic ring) on C atoms that are clearly in atomic rings. OpenBabel reports 37.
- The amount of vdw/ES interactions are inconsistent, Julia Ball has slightly more.


### Problem: atom types
- many atom types are being perceived wrong or not at all. see xxxxxxx.Txt for some examples.
- have started work on fixing it. The goal is to make it conform to Openbabel
- may use GetType in forcefieldmmff94.cpp ll 1481-2600
- solution to some of it: double bonds in aromatic rings arent double bonds "=", they are aro bonds ":"
- a false double bond between atoms num 7 and 8
- solved

#### BALL program workflow
function test(name)
                  syst = load_sdfile(ball_data_path("../test/data/" * name))        
                  mm= MMFF94FF(syst)
                  for mol in molecules(syst)

                  foreach(println, parse.(Int, (atoms_df(mol).atom_type)))

                  println("\n\n")
                  end
                  end

export BABEL_DATADIR=/usr/local/share/openbabel
cd /mnt/c/users/Dan/babelprogs

cd /mnt/c/users/Dan/ballprogs/build
export BALL_DATA_PATH=/mnt/c/Users/Dan/wslball/ball/data
cd build
cmake ..
make
./BALLprogs

# Links
[MMFF 94 paper](https://onlinelibrary.wiley.com/doi/10.1002/(SICI)1096-987X(199604)17:5/6%3C490::AID-JCC1%3E3.0.CO;2-P)

G:\Downloads\J Comput Chem - April 1996 - Halgren - Merck molecular force field  I  Basis  form  scope  parameterization  and (1)

BALL docu - file:///G:/Program%20Files/BALL-1.4.79/share/BALL/doc/BALL/index.html

[BALL git](https://github.com/BALL-Project/ball)

[Ball.jl git](https://github.com/hildebrandtlab/BiochemicalAlgorithms.jl/tree/ah/amberff)

[Ball.jl docu](./docs/build/index.html)

BALL Folder - G:\Program Files (x86)\ball\source

MMFF94.C Line 140 starts init

# Master Notes

## Forcefield  
- MMFF94: Read paper to find out materials, procedure

## MMFF94
- pairwise function that considers every pair
- 7 terms that are summed up
    - Bend
    - StretchBend
    - Stretch
    - VDW
    - ES
    - Tors
    - OOP



# Master - Outline

## Intro
