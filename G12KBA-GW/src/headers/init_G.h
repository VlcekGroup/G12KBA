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

#ifndef init_G_h
#define init_G_h
#include "HF.h"
#include "mkl_lapacke.h"
#include <vector>
#include <complex>
#include <algorithm>
#include <iostream>
using namespace std;

void eigenstates(vector<complex<double>> &h1, vector<complex<double>> &evecs,vector<float> &evals, int Ns,int Nb)
{
	MKL_INT n = Ns*Nb, info;
    vector<MKL_Complex8> a(h1.size());
    
    for (int j=0; j<n; ++j){
        for(int i = j; i<n;++i){
            a[i*n + j] = {h1[i*n+j].real(),h1[i*n + j].imag()};
        }
    }
    info = LAPACKE_cheev( LAPACK_ROW_MAJOR, 'V', 'L', n,& *a.begin(), n, & *evals.begin());
    for (int j=0; j < n; ++j){
        for(int i = 0; i < n; ++i){
            evecs[i*n + j] = {a[i*n+j].real,a[i*n + j].imag};
        }
    }
//    float temp1 = 0;
//    float temp2 = 0;
//    for(int i = 0; i<int(Ns/2);++i){
//    	temp1+=2*evals[i];
//	
//    }
//    temp2 = temp1-evals[int(Ns/2)-1] + evals[int(Ns/2)]; 
//    cout<<temp1<<" "<<temp2<<" "<<abs(temp1-temp2)<<"\n";
    if(info!=0)
        cout<<"Error in diagonalization\n";
}


void make_G0(vector<complex<double>> &G1u, int Ns, int Nb,int p_number,double U, double t1, double t2, double t3, vector<double> &epsilon)
{
    vector<complex<double>> hu(Ns*Ns*Nb*Nb);
    vector<double> V;
    h_HFx_quench(hu,G1u, G1u, Ns, Nb, U,t1,t2,0,epsilon,V,0);
	complex<double> im = {0,1};
	vector<complex<double>> evecs_u(Ns*Nb*Ns*Nb);
    vector<float> evals_u(Ns*Nb);
    eigenstates(hu,evecs_u,evals_u,Ns,Nb);

    for(int i = 0; i<Ns*Nb; ++i){
	for(int j = 0; j < Ns*Nb; ++j){
            G1u[i*Ns*Nb + j] = 0;
	    for(int k = 0;k<p_number;++k){
		G1u[i*Ns*Nb + j] += evecs_u[i*Ns*Nb + k]*conj(evecs_u[j*Ns*Nb+k])*im;
            }
	}
    }
}


void init_G_test(vector<complex<double>> &G1u, int Ns, int Nb)
{
    for(int i = 0; i < int(Ns/2); ++i){
        for(int a = 0; a < int(Nb); ++a){
            G1u[(i+a*Ns)*Nb*Ns + i + a*Ns] = im;
        }
    }
}
#endif





