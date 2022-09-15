#include "kron_delta.h"
#include <vector>
#include <complex>
using namespace std;

void h_HF_no_pulse(vector<complex<double>> &h, vector<complex<double>> &G1, int Ns, int Nb, double U, double t1, double t2,vector<double> epsilon)
{
        //Function that updates the Hartree-Fock Hamiltonian
	complex<double> im = {0,1};
    for(int i = 0; i < Ns; ++i){
        for(int a = 0; a < Nb; ++ a){
            h[(i + a*Ns)*Nb*Ns + (i + a*Ns)] = epsilon[a];
            h[((i+1)%Ns + a*Ns)*Nb*Ns + (i + a*Ns)] = -t1;
            h[(i + a*Ns)*Nb*Ns + ((i+1)%Ns + a*Ns)] = -t1;
            if(a < Nb-1){
                h[(i+a*Ns)*Ns*Nb + (i + (a+1)*Ns)] = - t2;
                h[(i+(a+1)*Ns)*Ns*Nb + (i + a*Ns)] = - t2;
            }
            for(int g = 0; g < Nb; ++g){
                if(g!=a){
                    h[(i +a*Ns)*Nb*Ns + (i+a*Ns)] += -im*U*G1[(i+g*Ns)*Ns*Nb + i + g*Ns];
                }
            }
            for(int b = 0; b < Nb; ++b){
                h[(i +a*Ns)*Nb*Ns + (i+b*Ns)] += im*U*kron_delta_bar(a,b)*G1[(i+b*Ns)*Ns*Nb + i + a*Ns];
            }
        }
    }
}

void h_HF_wavelet_pulse(vector<complex<double>> &h, vector<complex<double>> &G1,double t, int Ns, int Nb, double U, double t1, double t2,double P, double sig, double w,vector<double> epsilon)
{
        //Function that updates the Hartree-Fock Hamiltonian and has the inclusion of a Gaussian wavelet
	complex<double> im = {0,1};
    for(int i = 0; i < Ns; ++i){
        for(int a = 0; a < Nb; ++ a){
            h[(i + a*Ns)*Nb*Ns + (i + a*Ns)] = epsilon[a];
            h[((i+1)%Ns + a*Ns)*Nb*Ns + (i + a*Ns)] = -t1;
            h[(i + a*Ns)*Nb*Ns + ((i+1)%Ns + a*Ns)] = -t1;

            double eps = .01;
            if(a < Nb-1){
                h[(i+a*Ns)*Ns*Nb + (i + (a+1)*Ns)] = P*exp(-pow(t/eps,2))*sin(w*t) - t2;
                h[(i+(a+1)*Ns)*Ns*Nb + (i + a*Ns)] = P*exp(-pow(t/eps,2))*sin(w*t) - t2;
            }
            for(int g = 0; g < Nb; ++g){
                if(g!=a){
                    h[(i +a*Ns)*Nb*Ns + (i+a*Ns)] += -im*U*G1[(i+g*Ns)*Ns*Nb + i + g*Ns];
                }
            }
            for(int b = 0; b < Nb; ++b){
                h[(i +a*Ns)*Nb*Ns + (i+b*Ns)] += im*U*kron_delta_bar(a,b)*G1[(i+b*Ns)*Ns*Nb + i + a*Ns];
            }
        }
    }

}

void h_HF_delta_kick(vector<complex<double>> &h, vector<complex<double>> &G1, int Ns, int Nb, double U, double t1, double t2,double P, double sig, double w,vector<double> epsilon)
{
	//Still needs to be completed
        //Function that updates the Hartree-Fock Hamiltonian and has the inclusion of a Gaussian wavelet  
	complex<double> im = {0,1};

    for(int i = 0; i < Ns; ++i){
        for(int a = 0; a < Nb; ++ a){
            h[(i + a*Ns)*Nb*Ns + (i + a*Ns)] = epsilon[a];
            h[((i+1)%Ns + a*Ns)*Nb*Ns + (i + a*Ns)] = -t1;
            h[(i + a*Ns)*Nb*Ns + ((i+1)%Ns + a*Ns)] = -t1;

            if(a < Nb-1){
                    h[(i+a*Ns)*Ns*Nb + (i + (a+1)*Ns)] = -t2;
                    h[(i+(a+1)*Ns)*Ns*Nb + (i + a*Ns)] = -t2;
            }
            for(int g = 0; g < Nb; ++g){
                if(g!=a){
                    h[(i +a*Ns)*Nb*Ns + (i+a*Ns)] += -im*U*G1[(i+g*Ns)*Ns*Nb + i + g*Ns];
                }
            }
            for(int b = 0; b < Nb; ++b){
                h[(i +a*Ns)*Nb*Ns + (i+b*Ns)] += im*U*kron_delta_bar(a,b)*G1[(i+b*Ns)*Ns*Nb + i + a*Ns];
            }
        }
    }
}




