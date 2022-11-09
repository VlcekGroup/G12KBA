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


#ifndef collision_integral_h
#define collision_integral_h
#include "kron_delta.h"
#include <vector>
#include <complex>
using namespace std;


void collision_int(vector<complex<double>> &Ix, vector<complex<double>> &G2xyxy, vector<complex<double>> &G2xxxx, int Ns, int Nb,double U, vector<double> &V)
{
    //Function to calculate the collision integral from the 2 particle greens function.  Implements equation 6 in corresponding note/equation 31 in Bonitz paper
    int Ns2 = int(pow(Ns,2));
        int Nb2 = int(pow(Nb,2));
        complex<double> im = {0,1};
#define IDX_4D(x,y,z,t)  ((x) + (y)*Ns + ((z) + (t)*Ns)*Nb*Ns)
#define IDX_8D(x,y,z,t,m,n,q,w) ((x) + (y)*Ns + ((z) + (t)*Nb)*Ns2 + ((m)+(n)*Ns)*Ns2*Nb2 + ((q) + (w)*Nb)*Ns2*Ns2*Nb2)
        
    for(int i = 0; i<Ns; ++i){
        for(int j = 0; j<Ns; ++j){
            for(int a = 0; a<Nb; ++a){
                for(int b = 0; b<Nb; ++b){
                    int idx1 = IDX_4D(j,b,i,a);
                    int idx2 = IDX_8D(j,i,b,a,i,i,a,a);
                    Ix[idx1] = -im*U*G2xyxy[idx2];
                    for(int m = 0; m<Nb;++m){
                        if(m != a){
                            int idx3 = IDX_8D(j,i,b,m,i,i,a,m);
                            Ix[idx1] += -im*U*G2xxxx[idx3];
                            Ix[idx1] += -im*U*G2xyxy[idx3];
                        }
                    }
                    
                    for(int n = 1; n < V.size()+1; ++n){
                        if(i+n < Ns){
                            for(int g = 0; g < Nb; ++g){
	                            int idx4 = IDX_8D(j,i+n,b,g,i,i+n,a,g);
                                Ix[idx1] += -im*V[n-1]*(G2xyxy[idx4] + G2xxxx[idx4]);
                            }
                        }
                        if(i-n >= 0){
                            for(int g = 0; g < Nb; ++g){
                            int idx5 = IDX_8D(j,i-n,b,g,i,i-n,a,g);
                                Ix[idx1] += -im*V[n-1]*(G2xyxy[idx5] + G2xxxx[idx5]);
                            }
                        }
                    }
                }
            }
        }
    }
}
#endif

