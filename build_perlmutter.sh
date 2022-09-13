module unload darshan
module load PrgEnv-nvidia cudatoolkit
module load cray-libsci
cd src
make clean
make G1_G2_main

