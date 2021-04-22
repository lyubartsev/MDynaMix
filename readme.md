#       MDynaMix: a molecular dynamics program                 
##               v.5.3.0     

===================================================================

## Authors:

Alexander Lyubartsev     (alexander.lyubartsev@mmk.su.se)
Aatto Laaksonen          (aatto@mmk.su.se)
Alexei Nikitin           (agilemolecule@gmail.com)  

Div. Physical Chemistry
Department of Materials and Environmental Chemistry
Stockholm University
S 106 91   Stockholm
Sweden

=====================================================================

## Description:

MDynaMix is a general purpose molecular dynamics code for simulations
of mixtures  of either rigid or flexible molecules, interacting by
AMBER-like force field in a periodic rectangular cell. In the case of 
flexible molecules the double time step algorithm is used.
Algotithms for NVE, NVT, NPT and anisotropic NPT are implemented, 
Ewald sum for treatment of the electrostatic interactions,
account for quantum effects using path integral approach,
calculation of solvation free energy by expanded ensemble method. 
The program can be run both in sequential and parallel execution, in
the latter case the MPI parallel environment is required.    

=======================================================================

## Quick compilation:

$ make

or

$ make mpi


See futher details in "doc" directory, or upload repository mdynamix-manual

=======================================================================