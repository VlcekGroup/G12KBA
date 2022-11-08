# G12KBA
HF-GKBA code with efficient time-linear formulation employing one- and two-body Green's functions

# Overview:
Efficient and scalable c++ code, optimized for HW operated by the National Energy Research Scientific Computing Center, implements the propagation of the equal-time non-equilibrium one- and two-body Green's function (NEGF).  Part of the methodology uses the "G1-G2" scheme presented in https://journals.aps.org/prb/abstract/10.1103/PhysRevB.101.245101. The code prepares the stationary state via adiabatic switching procedure, employs the GW self-energy, and propagates the equal time NEGF for lattice like models with multiple orbitals per site and long range interactions.  

This implementation allows for correlated initial states to be prepared from the non-interacting ground state of the Hamiltonian; the driving from equilibrium is carried out via quenching (rapidly changing the onsite energy of certain sites), other out of equilibrium preperations can be implemented by editing HF.h. 


# Required packages:
MKL Lapack is used for matrix operations however cblas can easily be subsituted

Optimized with openMP for parallel execution
# Compilation and Execution:

Source mkl libraries ie "source /path/to/mkl.sh"

Compile with "g++ <filename>.cpp -fopenmp -lmkl_intel_lp64 -lmkl_gnu_thread -lpthread -lm -ldl -lmkl_core -O3 -o <output_file>"
or run make command
  
Run with ./<output_file>
 
Alternatively run cmd ./run.sh
 
# Output:
Code outputs the full one particle Green's function along the entire calculated trajectory as well as the dipole to .txt file.  the two particle Green's function is also readily produced from the code.

# Additional Notes
## Future Versions: 

- Reconstruction of off diagonal components of Green's function
- More general Hamiltonians(different driving fields)
- Functions for calculating more quantities from G1 and G2(spectrum, energies, etc)

 # Acknowledgments:
 
 This material is based upon work supported by the U.S. Department of Energy, Office of Science, Office of Advanced Scientific Computing Research and Office of Basic Energy Sciences, Scientific Discovery through Advanced Computing (SciDAC) program under Award Number DE-SC0022198
