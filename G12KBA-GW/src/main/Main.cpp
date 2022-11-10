
//
//
//   The routine(s) in this file are a part of the
//                     G12KBA
//   suite, developed 2022, and copyrighted
//   to the authors: Cian Reeves and Vojtech Vlcek
//   at the University of California, Santa Barbara
//   and Khaled Ibrahim
//   at Lawrence Berkeley National Lab, Berkeley.
//
//
//  If you use or modify any part of this routine
//  the header should be kept and unmodified.
//
//
//

#include "../headers/write_to_file.h"
#include "../headers/RK_solver.h"
#include "../headers/AS.h"
#include "../headers/init_G.h"
#include "../headers/dipole.h"
#include "../headers/read_input.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <complex>
#include <chrono>
#include <string>
#include <ctime>
using namespace std;
complex<double> im={0,1};
extern int Ns;
extern int Nb;
extern int Ns2;
extern int Nb2;
extern double t1;
extern double t2;
extern double t3;
extern vector<double> epsilon;

int main()
{
    
    assign_vals();
    vector<complex<double>> G1x(Ns*Ns*Nb*Nb);
    vector<complex<double>> G2xxxx(Ns2*Ns2*Nb2*Nb2);
    vector<complex<double>> G2xydy(Ns2*Ns2*Nb2*Nb2);
    vector<vector<complex<double>>> G1x_diag;
    make_G0(G1x, Ns, Nb,int(Ns/2),0, t1, t2, t3, epsilon);
    Adiabatic_switching(G1x,G2xxxx,G2xydy,G1x_diag);
    RK4(G1x,G2xxxx,G2xydy,G1x_diag);
    vector<double> p(G1x_diag.size());
    dipole(p,G1x_diag);
    write_to_cfile2D("G1.txt",G1x_diag);
    write_to_rfile("dipole.txt",p);

    
    return 0;
}

