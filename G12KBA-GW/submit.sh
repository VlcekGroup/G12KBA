#!/bin/bash
home=`pwd`
cd src
source /opt/intel/parallel_studio_xe_2019.0.045/psxevars.sh
export OMP_NUM_THREADS=4
Ns=4
quench_type=full
EHM=true
for i in 0.50
do
        for U in 0.1 0.3 0.5 1.0 1.5
        do
		if [ "$EHM" == "true" ]; 
		then

                	mkdir "${quench_type}_quench_EHM_Ns_${Ns}_U_${U}_extent_${i}"; cd "${quench_type}_quench_EHM_Ns_${Ns}_U_${U}_extent_${i}"
		else		
	        	mkdir "${quench_type}_quench_onsite_Ns_${Ns}_U_${U}_extent_${i}"; cd "${quench_type}_quench_onsite_Ns_${Ns}_U_${U}_extent_${i}"
		fi
                cp $home"/src/main/Main.cpp" .
                cp $home"/src/main/Makefile" .

                sed -i "s/double\ interaction_strength=.*/double\ interaction_strength=$U;/" Main.cpp
                sed -i "s/int\ quench_extent=.*/int\ quench_extent=int(Ns*$i);/" Main.cpp
                sed -i "s/int\ Ns=.*/int\ Ns=$Ns;/" Main.cpp
                sed -i "s/bool\ Extended=.*/bool\ Extended=$EHM;/" Main.cpp
                sed -i "s/string\ quench_type=.*/string\ quench_type=\"$quench_type\";/" Main.cpp

                rm G1_G2_main
                make

                ./G1_G2_main >& log_$U &
                cd $home"/src/"
        done
done
cd ..
