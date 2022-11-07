# G12KBA
G1-G2 scheme for the HF-GKBA

# Overview:
c++ code which implements the G1-G2 scheme for the propagation of of the equal time non-equilibrium Green's function (NEGF).  The method implemented is based off the following paper https://journals.aps.org/prb/abstract/10.1103/PhysRevB.101.245101. This code uses the GW self-energy and allows to the propagation of the equal time NEGF for Hubbard like models with multiple orbitals per site and long range interactions.  

This implementation allows for correlated initial states to be prepared from the non-interacting ground state of the Hamiltonian through Adiabtatic switching. Currently driving from equilibrium is carried out via quenching(rapidly changing the onsite energy of certain sites), other out of equilibrium preperations can be implemented by editing HF.h. 


# Required packages:
MKL Lapack is currently used for matrix operations however cblas can easily be subsituted
The code is optimized for parallel execution using openMP

# Compilation and Execution:

Source mkl libraries e.g "source /opt/intel/parallel_studio_xe_2019/psxevars.sh"

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
