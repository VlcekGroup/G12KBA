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
    float temp1 = 0;
    float temp2 = 0;
    for(int i = 0; i<int(Ns/2);++i){
    	temp1+=2*evals[i];
	
    }
    temp2 = temp1-evals[int(Ns/2)-1] + evals[int(Ns/2)]; 
    cout<<temp1<<" "<<temp2<<" "<<abs(temp1-temp2)<<"\n";
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


void init_G(vector<complex<double>> &G1u, int Ns, int Nb)
{
    for(int i = 0; i < int(Ns/2); ++i){
        for(int a = 0; a < int(Nb); ++a){
            G1u[(i+a*Ns)*Nb*Ns + i + a*Ns] = im;
        }
    }
}



void init_G_eq(vector<complex<double>> &G1u, int Ns, int Nb)
{
	for(int i = 0; i < int(Ns); ++i){
		for(int j = 0; j <  Ns; ++j){
 	       		for(int a = 0; a < int(Nb); ++a){
				G1u[(i+a*Ns)*Nb*Ns + j + a*Ns] = .5*im;
			}
        	}
    	}
}


void init_G_Ns_10(vector<complex<double>> &G1u,vector<complex<double>> &G1d, int Ns, int Nb)
{
vector<vector<complex<double>>> temp = {
{ {0,0.5},{0,0.428473},{0,-1.38778e-16},{0,-0.178177},{0,-3.747e-16},{0,0.122529},{0,2.08167e-17},{0,-0.100305},{0,9.71445e-17},{0,0.0916312} },
{ {0,0.428473},{0,0.5},{0,0.250296},{0,-4.71845e-16},{0,-0.0556478},{0,-2.77556e-17},{0,0.0222236},{0,2.498e-16},{0,-0.00867433},{0,3.46945e-17} },
{ { 0,-1.11022e-16},{0,0.250296},{0,0.5},{0,0.372825},{0,-3.60822e-16},{0,-0.155953},{0,-7.80626e-17},{0,0.113855},{0,6.93889e-17},{0,-0.100305} },
{ {0,-0.178177},{0,-4.71845e-16},{0,0.372825},{0,0.5},{0,0.27252},{0,-2.08167e-16},{0,-0.0643221},{0,-3.98986e-16},{0,0.0222236},{0,1.04083e-16} },
{ {0,-4.02456e-16},{0,-0.0556478},{0,-3.60822e-16},{0,0.27252},{0,0.5},{0,0.364151},{0,-2.498e-16},{0,-0.155953},{0,-3.40006e-16},{0,0.122529} },
{ {0,0.122529},{0,-4.16334e-17},{0,-0.155953},{0,-2.08167e-16},{0,0.364151},{0,0.5},{0,0.27252},{0,0},{0,-0.0556478},{0,-4.16334e-16} },
{ {0,3.46945e-17},{0,0.0222236},{0,-5.72459e-17},{0,-0.0643221},{0,-2.498e-16},{0,0.27252},{0,0.5},{0,0.372825},{0,8.32667e-17},{0,-0.178177} },
{ {0,-0.100305},{0,2.77556e-16},{0,0.113855},{0,-3.7817e-16},{0,-0.155953},{0,0},{0,0.372825},{0,0.5},{0,0.250296},{0,2.498e-16} },
{ {0,8.32667e-17},{0,-0.00867433},{0,6.93889e-17},{0,0.0222236},{0,-3.40006e-16},{0,-0.0556478},{0,5.55112e-17},{0,0.250296},{0,0.5},{0,0.428473} },
{ {0,0.0916312},{0,4.85723e-17},{0,-0.100305},{0,9.02056e-17},{0,0.122529},{0,-4.02456e-16},{0,-0.178177},{0,2.22045e-16},{0,0.428473},{0,0.5} }
};

for(int i = 0; i < int(Ns); ++i){
                for(int j = 0; j <  Ns; ++j){
                        for(int a = 0; a < int(Nb); ++a){
                                G1u[(i+a*Ns)*Nb*Ns + j + a*Ns] = temp[i][j];
                                G1d[(i+a*Ns)*Nb*Ns + j + a*Ns] = temp[i][j];
                        }
                }
        }
}

void init_G_Ns_8(vector<complex<double>> &G1u, int Ns, int Nb){
	vector<vector<complex<double>>> temp ={ 
{{0,0.5},{0,0.430899},{0,-8.358e-17},{0,-0.183363},{0,1.10812e-16},{0,0.131279},{0,3.88977e-17},{0,-0.114262}},
{{0,0.430899},{0,0.5},{0,0.247536},{0,1.76007e-17},{0,-0.0520841},{0,-4.98684e-17},{0,0.0170171},{0,-1.76676e-17}},
{{0,-1.25213e-16},{0,0.247536},{0,0.5},{0,0.378815},{0,-4.28942e-17},{0,-0.166346},{0,-1.19352e-16},{0,0.131279}},
{{0,-0.183363},{0,1.76007e-17},{0,0.378815},{0,0.5},{0,0.264553},{0,6.19115e-17},{0,-0.0520841},{0,4.82064e-18}},
{{0,1.10812e-16},{0,-0.0520841},{0,-4.28942e-17},{0,0.264553},{0,0.5},{0,0.378815},{0,8.11207e-17},{0,-0.183363}},
{{0,0.131279},{0,-2.21128e-17},{0,-0.166346},{0,6.19115e-17},{0,0.378815},{0,0.5},{0,0.247536},{0,-2.16475e-17}},
{{0,4.58366e-17},{0,0.0170171},{0,-1.60986e-16},{0,-0.0520841},{0,8.11207e-17},{0,0.247536},{0,0.5},{0,0.430899}},
{{0,-0.114262},{0,-1.07287e-17},{0,0.131279},{0,4.82064e-18},{0,-0.183363},{0,-7.76968e-18},{0,0.430899},{0,0.5}}
};
for(int i = 0; i < Ns; ++i){
   for(int j = 0; j <  Ns; ++j){
      for(int a = 0; a < Nb; ++a){
                                G1u[(i+a*Ns)*Nb*Ns + j + a*Ns] = temp[i][j];
                        }
                }
        }
}
