#!/bin/bash
cd src/main
rm -f G1_G2_main
export OMP_NUM_THREADS=4
source /opt/intel/parallel_studio_xe_2019.0.045/psxevars.sh 
make
./G1_G2_main
cd ..
