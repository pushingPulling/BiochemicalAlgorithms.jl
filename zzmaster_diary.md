## MMFF94 description

in MMFF94 paper sind die params und die Energy func beschrieben. [paper](https://onlinelibrary.wiley.com/doi/10.1002/(SICI)1096-987X(199604)17:5/6%3C490::AID-JCC1%3E3.0.CO;2-P)

Es sind die atom types gelistet

MMFF94 ist für MDS geeignet, nicht für energy min

## Forcefields in Biochem

`AbstractForceFieldComponent` is for BondAngles, BondBends, torsion etc
`ForceField` enthält ein system, parameter und mehrere `AbstractForceFieldComponent`
in amber.jl wird ein FF erstellt. Das sollte das Ziel sein
Es gibt `compute_forces` und `compute_energy` was man auf die FF einwirkt.

# Weiteres Vorgehen

#### Finde beispiel bei dem die Parameter zuweisung gehen. analysiere dann, was genau in dem FF gespeichert ist.
 
Kann mir ah_amberff anschauen um zu gucken wie ich MMFF machen würde
#### Versuche einen Ablauf an funktionen zu finden (e.g. FF -> Components zu FF hinzufügen -> compute forces)
-   alle components durchschauen
-   atomtype verstehen
-   forcefield.jl durcharbeiten
-   stuff in src/forcefields/AMBER durcharbeiten
BALL implement von MMF94 nachbauen

# Implement diary
-   füge molgraph 0.14.2 hinzu wegen kekulization
    -   interface hat sich geändert, graphs werden jetzt benutzt. hab imprinzip nur am umschreiben (aka einarbeiten) gearbeitet.
-   abgebrochen, zu inkonsistent und das ganze interface hat sich bei Molgraph 0.13 -> 0.14 geändert
- mache ab jetzt mit molgraph@0.13 weiter
- implementiere c++ BALL'S kekulizer, aromatizer
- wtf macht l179 aromaticiyProcessor
    - debug/look in c++ BALL?
#### BALL program workflow

in /mnt/c/users/Dan/ballprogs
export BALL_DATA_PATH=/mnt/c/Users/Dan/wslball/ball/data
cd build
cmake ..
make
./BALLprogs

# Links
[MMFF 94 paper](https://onlinelibrary.wiley.com/doi/10.1002/(SICI)1096-987X(199604)17:5/6%3C490::AID-JCC1%3E3.0.CO;2-P)

BALL docu - file:///G:/Program%20Files/BALL-1.4.79/share/BALL/doc/BALL/index.html

[BALL git](https://github.com/BALL-Project/ball)

[Ball.jl git](https://github.com/hildebrandtlab/BiochemicalAlgorithms.jl/tree/ah/amberff)

[Ball.jl docu](./docs/build/index.html)

BALL Folder - G:\Program Files (x86)\ball\source

MMFF94.C Line 140 starts init


## Gitarre
Master of Puppies YEAAAAAAAAAAAA
Tamacun
holy war
rasputin
jojo golden wind