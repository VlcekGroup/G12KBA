# G12KBA
HF-GKBA code with efficient time-linear formulation employing one- and two-body Green's functions

# Overview:
Efficient and scalable c++ code, optimized for HW operated by the National Energy Research Scientific Computing Center, implements the propagation of the equal-time non-equilibrium one- and two-body Green's function (NEGF).  Part of the methodology uses the "G1-G2" scheme presented in https://journals.aps.org/prb/abstract/10.1103/PhysRevB.101.245101. The code prepares the stationary state via adiabatic switching procedure, employs the GW self-energy, and propagates the equal time NEGF for lattice like models with multiple orbitals per site and long range interactions.  

This implementation allows for correlated initial states to be prepared from the non-interacting ground state of the Hamiltonian; the driving from equilibrium is carried out via quenching (rapidly changing the onsite energy of certain sites), other out of equilibrium preperations can be implemented by editing HF.h. 


# Required packages & files:
MKL Lapack is used for matrix operations however cblas can easily be subsituted

Optimized with openMP for parallel execution

input file called "input" must be included in same directory as main.cpp.  If not included default values will be assigned to variables.  sample input is included.
# Compilation and Execution:

Source mkl libraries ie "source /path/to/mkl.sh"

Compile with "g++ <filename>.cpp -fopenmp -lmkl_intel_lp64 -lmkl_gnu_thread -lpthread -lm -ldl -lmkl_core -O3 -o <output_file>"
or run make command
  
Run with ./<output_file>
 
Alternatively run cmd ./run.sh
 
# Output:
Code outputs the full one particle Green's function along the entire calculated trajectory as well as the dipole to .txt file.  the two particle Green's function is also readily produced from the code.


# Additional Notes
## Input Parameters:
- HF = true/false, determines whether HF or full HF-GKBA calculation will be performed(default=false)
- EHM true/false, determines if there are long range interactions or onsite interactions only
- t1 = real number, intersite hopping strength
- t2 = real number, interband hopping strength
- U = real number, onsite two-body interaction strength
- time_steps = positive integer, number of time steps to be taken by RK stepper
- dt = positive real number, step size for RK solver
- quench_strength = real number, magnitude of system quench
- Ns = positive integer, number of sites in the chain
- Nq = positive integer, number of contiguous sites to be quenched counting from the first
- Nb = positive integer, number of orbitals per site
- q_type = "full","pulse","none", type of quench on system.  full lasts for all time, quench has finite time width
- AS_rate = positive real number, rate at which interaction is turned on
- AS_midpoint = positive real number, time at which AS function = .5
- quench_rate = positive real number, rate at which quench is turned on
- quench_on_midpoint = positive real number, time at which quench .5*quench_strength and is  increasing
- quench_off_midpoint = positive real number, time at which quench .5*quench_strength and is  decreasing
## Spin symmetry:
Spin symmetry is assumed so only one spin species is time evolved

## Future Versions: 

- Reconstruction of off diagonal components of Green's function
- More general Hamiltonians(different driving fields)
- Functions for calculating more quantities from G1 and G2(spectrum, energies, etc)

 # Acknowledgments:
 
 This material is based upon work supported by the U.S. Department of Energy, Office of Science, Office of Advanced Scientific Computing Research and Office of Basic Energy Sciences, Scientific Discovery through Advanced Computing (SciDAC) program under Award Number DE-SC0022198
