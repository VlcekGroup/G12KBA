G1-G2-Scheme
G1-G2 scheme for the HF-GKBA

cpp code which implements the G1-G2 scheme for the Hubbard type hamiltonian with multiple orbitals/spins per site.  Interband couplings are allowed in this implementation. The lack of interband/interspin transiton in the standard Hubbard model allows
for block diagonal structure to be taken advantage of.  In this code full NbxNs and (NbxNs)^2 matrices are used  where Nb is the number of bands/spins per site and Ns is the number of sites in the system

The G1-G2 scheme is an approximate method for the time propagation of non equilibrium Green's functions.  The method implemented is based off the following paper https://journals.aps.org/prb/abstract/10.1103/PhysRevB.101.245101.

The code allows for the inclusion of a driving field that allows for transitions between different orbitals. Currently only the code only supports propagation along the time diagonal. This code uses the GW self energy

Future implementations will include 
-Reconstruction of off diagonal components of Green's function
-More general Hamiltonians(Longer range couplings, different driving fields)

Compilation:

Source mkl libraries e.g "source /opt/intel/parallel_studio_xe_2019/psxevars.sh"

Compile with "g++ <filename>.cpp -fopenmp -lmkl_intel_lp64 -lmkl_gnu_thread -lpthread -lm -ldl -lmkl_core -O3 -o <output_file>"

  
 Run with ./<output_file>
