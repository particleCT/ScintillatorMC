# ScintillatorMC

Project to mimic a scintilaltor based detector for proton radiography.

Compile by running source compile.sh

Prerequesite:
- Geant4
- ROOT

An example command to run is
./bin/Linux-g++/ScintillatorMC NParticle Energy Model Angle Thickness Thread ANumber


where you replace every variable by what you wishes (e.g. replace Energy by 200).
Except if you are using an XCAT phantom, where you run with

./bin/Linux-g++/ScintillatorMC NParticle Energy XCAT Angle Thickness Thread ANumber LungPhantom/choose_a_phase.root