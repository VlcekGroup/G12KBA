# G1-G2-Scheme
G1-G2 scheme for the HF-GKBA

cpp code which implements the G1-G2 scheme for the Hubbard type hamiltonian with multiple orbitals/spins per site. The G1-G2 scheme is an approximate method for the time propagation of non equilibrium Green's functions.  The method implemented is based off the following paper https://journals.aps.org/prb/abstract/10.1103/PhysRevB.101.245101.

The code allows for the inclusion of a driving field that allows for transitions between different orbitals. Currently only the code only supports propagation along the time diagonal.  Only the second order Born self energy is currently implemented

Future implementations will include 
-Reconstruction of off diagonal components of Green's function
-GW self energies
-More general Hamiltonians(Longer range couplings, different driving fields)
