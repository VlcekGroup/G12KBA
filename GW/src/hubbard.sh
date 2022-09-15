#!/bin/sh
#SBATCH --nodes=8
#SBATCH --time=5
#SBATCH --constraint=knl
#SBATCH --qos=<QOS>
#SBATCH --account=m4022

export OMP_NUM_THREADS=8
srun -n <num_mpi_processes> -c <cpus_per_task> a.out

# perform any cleanup or short post-processing here
