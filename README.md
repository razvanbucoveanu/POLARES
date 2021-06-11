# POLARES: a C++ program for the simulation and event generation of elastic scattering of polarized leptons off nuclei

POLARES is a C++ program for elastic lepton-nucleon scattering with longitudinally polarized
electron beams. It includes QED radiative corrections at the one and two-loop level for
unpolarized, as well as for polarized incident leptons. It can be used as an integrator to
calculate cross sections and asymmetries for given kinematic conditions. It can also be
used as an event generator. The design of the program POLARES was developed in such a
way that it can be easily combined with the detector simulation software of the experiment.
It contains an option to generate events with varying energy of the incoming lepton. In
addition, we have included the possibility to switch between electron and muon as the
incident particle and between proton and carbon-12 as the target particle. Note, however,
that in the present version the two-loop corrections are based on an calculation which
assumes the lepton mass to be negligibly small compared with the momentum transfer.
This approximation may not be valid for muon scattering at low energies.

For more details about the program and instructions on how to install and use it please consult 
the manual found in the /doc folder. 

The present version is preliminary. 
