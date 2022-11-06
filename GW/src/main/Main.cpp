/*
*/
//User defined headers
#include "../headers/write_to_file.h"
#include "../headers/RK_solver.h"
#include "../headers/print_matrix.h"//Prints matrices (used for testing)
#include "../headers/init_G.h"
#include "../headers/dipole.h"
//#include <unistd.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <complex>
#include <chrono>
#include <string>
#include <ctime>
using namespace std;

int main()
{
    int Nb=1;
    int Ns=4;
    int time_steps=105000;
    double step_size=.01;
    double prop_time = double(time_steps)*step_size;
    double inter_site_hopping=1;
    double interaction_strength=1;
    double quench=1;
    double alpha=.7;
    bool Har_Fock=false;
    bool Extended=true;
    int quench_extent=int(Ns/4);
    string quench_type="pulse";
    hubbard_RK_solver GW(interaction_strength,inter_site_hopping,0,0,0,Ns,Nb,0,step_size,prop_time,quench,alpha,Extended,Har_Fock,quench_extent,quench_type);

    vector<complex<double>> G1u(GW.Ns*GW.Ns*GW.Nb*GW.Nb);
    vector<complex<double>> G2uuuu(GW.Ns2*GW.Ns2*GW.Nb2*GW.Nb2);
    vector<complex<double>> G2udud(GW.Ns2*GW.Ns2*GW.Nb2*GW.Nb2);
    vector<vector<complex<double>>> G1u_diag;
    make_G0(G1u, Ns, Nb,int(Ns/2),0, GW.t1, GW.t2, GW.t3, GW.epsilon);
    GW.Adiabatic_switching(G1u,G2uuuu,G2udud,G1u_diag);
    cout<<"Adiabatic Switching complete: Initiating Quench\n";
    auto start = chrono::steady_clock::now();
    int check=GW.RK4(G1u,G2uuuu,G2udud,G1u_diag);
    vector<double> p(G1u_diag.size());

    dipole(p,G1u_diag);
    auto end = chrono::steady_clock::now();
    if(check==0){
	cout<<"\n Execution Succesful, elapsed time:"<< chrono::duration<double>(end - start).count()<< " sec \n";
    }
    if(check==1){
	cout<<"\n Execution Failed, unphysical number density, elapsed time:"<< chrono::duration<double>(end - start).count()<< " sec \n";
    }
	
 
   
    string filename = "G1.txt";
    write_to_cfile2D(filename,G1u_diag);
    filename="dipole.txt";
    write_to_rfile(filename,p);

    
    return 0;
}

