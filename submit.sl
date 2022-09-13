#!/bin/bash
#SBATCH -A m2956
#SBATCH -C cpu
#SBATCH -q debug
#SBATCH -t 0:30:00
#SBATCH -N 1
#SBATCH --mail-user=kzibrahim@lbl.gov


module unload darshan
module load PrgEnv-nvidia cudatoolkit
module load cray-libsci

#salloc --nodes 1 --qos interactive --time 01:00:00 --constraint cpu   --account=m2956
cd src

APP=G1_G2_main
export OMP_NUM_THREADS=8
cores=$(($OMP_NUM_THREADS * 2))


srun -n 1 -c $cores ./$APP


