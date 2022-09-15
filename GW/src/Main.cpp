/*
*/

//User defined headers
#include "write_to_file.h"
#include "RK_solver.h"
#include "print_matrix.h"//Prints matrices (used for testing)

//#include <unistd.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <complex>
#include <chrono>
using namespace std;



int main()
{
	auto start = chrono::steady_clock::now();
	hubbard_RK_solver SB(1,1,0,0,1,2,2,0);

	vector<complex<double>> G1(SB.Ns*SB.Ns*SB.Nb*SB.Nb);
	vector<complex<double>> G2(SB.Ns2*SB.Ns2*SB.Nb2*SB.Nb2);
	vector<double> n1;
	vector<vector<complex<double>>> G_diag;
	double avg_wall_t_per_loop = 0;

	SB.init_G(G1);
	SB.RK_steps(n1,G1,G2,G_diag,avg_wall_t_per_loop);
	write_to_rfile("n1.txt", n1);
    
	auto end = chrono::steady_clock::now();
	cout<<"Time stepping complete, elapsed time:"<< chrono::duration<double>(end - start).count()<< " sec \n";

	return 0;
}

