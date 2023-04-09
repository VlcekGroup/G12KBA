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

#ifndef HF_h
#define HF_h
#include "kron_delta.h"
#include <vector>
#include <iostream>
#include <complex>
using namespace std;

// void h_HFx(vector<complex<double>> &h, vector<complex<double>> &G1x, vector<complex<double>> G1y, int Ns, int Nb, double U, double t1, double t2,vector<double> epsilon)
// {
//         //Function that updates the Hartree-Fock Hamiltonian
//     complex<double> im = {0,1};
//     for(int i = 0; i < Ns; ++i){
//         for(int a = 0; a < Nb; ++a){
//             h[(i + a*Ns)*Nb*Ns + (i + a*Ns)] = epsilon[a];
// 		if(i < Ns-1){
//                 h[((i+1)+a*Ns)*Nb*Ns + (i+a*Ns)] += -t1;
//                 h[(i+a*Ns)*Nb*Ns + ((i+1)+a*Ns)] += -t1;
//             }
//             if(a < Nb-1){
//                 h[(i+a*Ns)*Ns*Nb + (i + (a+1)*Ns)] = - t2;
//                 h[(i+(a+1)*Ns)*Ns*Nb + (i + a*Ns)] = - t2;
//             }
//             h[(i +a*Ns)*Nb*Ns + (i+a*Ns)] += -im*U*G1y[(i +a*Ns)*Nb*Ns + (i+a*Ns)];
            
//             for(int g = 0; g < Nb; ++g){
//                 if(g!=a){
//                     h[(i +a*Ns)*Nb*Ns + (i+a*Ns)] += -im*U*G1x[(i+g*Ns)*Ns*Nb + i + g*Ns];
//                     h[(i +a*Ns)*Nb*Ns + (i+a*Ns)] += -im*U*G1y[(i+g*Ns)*Ns*Nb + i + g*Ns];

//                 }
//             }
//             for(int b = 0; b < Nb; ++b){
//                 h[(i +a*Ns)*Nb*Ns + (i+b*Ns)] += im*U*kron_delta_bar(a,b)*G1x[(i+b*Ns)*Ns*Nb + i + a*Ns];
//             }
//         }
//     }
// }

void h_HFx_quench(vector<complex<double>> &h, vector<complex<double>> &G1x, vector<complex<double>> G1y, int Ns, int Nb, double U, double t1, double t2,double w, vector<double> &epsilon, vector<double> &V,int quench_extent)
{
        //Function that updates the Hartree-Fock Hamiltonian
    complex<double> im = {0,1};

    for(int i = 0; i < Ns; ++i){
        for(int a = 0; a < Nb; ++a){
//             h[((i)+a*Ns)*Nb*Ns + (i+a*Ns)] += 2*pow(-1,i);
            h[((i)+a*Ns)*Nb*Ns + (i+a*Ns)] += epsilon[a];
            if(i<quench_extent){
                h[((i)+a*Ns)*Nb*Ns + (i+a*Ns)] += w;
            }

            if(i < Ns-1){
                h[((i+1)+a*Ns)*Nb*Ns + (i+a*Ns)] += -t1;
                h[(i+a*Ns)*Nb*Ns + ((i+1)+a*Ns)] += -t1;
            }

            if(a < Nb-1){
                h[(i+a*Ns)*Ns*Nb + (i+(a+1)*Ns)] += - t2;
                h[(i+(a+1)*Ns)*Ns*Nb + (i+a*Ns)] += - t2;
            }

            h[(i+a*Ns)*Nb*Ns + (i+a*Ns)] += -im*U*G1x[(i+a*Ns)*Nb*Ns + (i+a*Ns)]; //Hartree term for onsite and same orbital

            for(int g = 0; g < Nb; ++g){
                h[(i+a*Ns)*Nb*Ns + (i+g*Ns)] += im*U*kron_delta_bar(a,g)*G1x[(i+g*Ns)*Ns*Nb + (i+a*Ns)]; //Fock for onsite different orbitals
                if(g!=a){
                    h[(i+a*Ns)*Nb*Ns + (i+a*Ns)] += -2.0*im*U*G1x[(i+g*Ns)*Ns*Nb + (i+g*Ns)];//factor of 2 for spin degeneracy Hartree term for oniste different orbitals
                }
            }
            for(int n = 1; n < V.size() + 1; ++n){
                if(i + n < Ns){
                    for(int d = 0; d < Nb; ++d){
                        h[(i+a*Ns)*Nb*Ns + (i+a*Ns)] += -2.0*im*V[n-1]*(G1x[((i+n)+d*Ns)*Nb*Ns + ((i+n)+d*Ns)]);//factor of 2 for spin degeneracy  Hartree term for sites i+j = n
                        h[((i+n)+a*Ns)*Nb*Ns + (i+d*Ns)] += im*V[n-1]*G1x[((i+n)+a*Ns)*Nb*Ns + (i+d*Ns)]; //Fock term  for sites seperated by n
                        h[(i+a*Ns)*Nb*Ns + ((i+n)+d*Ns)] += im*V[n-1]*G1x[(i+a*Ns)*Nb*Ns + ((i+n)+d*Ns)]; //Fock term for sites seperated by n
                    }
                }
                if(i - n >= 0){
                    for(int d = 0; d < Nb; ++d){
                        h[(i+a*Ns)*Nb*Ns + (i+a*Ns)] += -2.0*im*V[n-1]*(G1x[((i-n)+d*Ns)*Nb*Ns + ((i-n)+d*Ns)]);//factor of 2 for spin degeneracy  Hartree term for sites |i-j| = n
                    }
                }
            }
        }
    }
}


#endif
