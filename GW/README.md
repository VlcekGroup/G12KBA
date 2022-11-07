G1-G2-Scheme
G1-G2 scheme for the HF-GKBA

A cpp code that implements the G1-G2 scheme for the propagation of the one and two particle Green's function within the HF-GKBA.  This method implemented is based off the following paper https://journals.aps.org/prb/abstract/10.1103/PhysRevB.101.245101.  The current implementation is for Hubbard clusters with onsite or long range interactions(EHM) and  Nb bands/orbitals per site and Ns sites.  This implementation is for the GW self-energy only.  

The code also allows for the preperation of correlated initial states through by first preparing the ground state of the non interacting model and the application of adiabatic switching.  Once prepared in equilibrium the system can be driven from equilibrium with a pulse of some temporal width.  Modifications to the quench or driving dynamics can be added in the HF header file.

The output files of the code are the dipole moment of the system over the time period as well as the full Green's function.  The code can also readily provide the two particle Green's function.  The code also allows for the inclusion of a driving field that allows for transitions between different orbitals. Currently only the code only calculates the time diagonal elements of the two-time Green's function.

Required libraries:
Matrix routines are done though MKL lapack although these can be replaced by cblas
Code is parallelized with openmp


Future implementations will include 
-Reconstruction of off diagonal components of Green's function
-More general Hamiltonian input
-More general driving procedures
-More ouput possibilities(Energy,spectrum, etc)


Compilation:

Source mkl libraries e.g "source /opt/intel/parallel_studio_xe_2019/psxevars.sh"

Compile with "g++ <filename>.cpp -fopenmp -lmkl_intel_lp64 -lmkl_gnu_thread -lpthread -lm -ldl -lmkl_core -O3 -o <output_file>"(Or run make)
Run with ./<output_file>
Or run
./run.sh



