//
//  pi.h
//  
//
//  Created by Cian Reeves on 9/13/22.
//

#ifndef pi_h
#define pi_h
#include <vector>
#include <complex>
#include "kron_delta.h"
using namespace std;

void Pi_xyxy(vector<complex<double>> &pi, vector<complex<double>> &G1y, vector<complex<double>> &G2xyxy,vector<complex<double>> &G2xxxx, double U, int Ns, int Nb,vector<double> &V)
{
    int Ns2 = int(pow(Ns,2));
    int Nb2 = int(pow(Nb,2));
    complex<double> im = {0,1};
#define IDX_4D(x,y,z,t)  ((x) + (y)*Ns + ((z) + (t)*Ns)*Nb*Ns)
#define IDX_8D(x,y,z,t,m,n,q,w) ((x) + (y)*Ns + ((z) + (t)*Nb)*Ns2 + ((m)+(n)*Ns)*Ns2*Nb2 + ((q) + (w)*Nb)*Ns2*Ns2*Nb2)
#if NOTHREADS
    #pragma omp parallel for collapse(8)
#endif
    for(int i = 0; i < Ns; ++i){
        for(int j = 0; j < Ns; ++j){
            for(int k = 0; k < Ns; ++k){
                for(int l = 0; l < Ns; ++l){
                    for(int a = 0; a < Nb; ++a){
                        for(int b = 0; b < Nb; ++b){
                            for(int g = 0; g < Nb; ++g){
                                for(int d = 0; d < Nb; ++d){
                                    int idx1 = IDX_8D(k,l,g,d,i,j,a,b);
                                    int idx2 = IDX_4D(l,d,j,b);
                                    pi[idx1] = 0;
                                    for(int n = 0; n < V.size() + 1; ++n){
                                        if(n==0){
                                            for(int r = 0; r < Nb; ++r){
                                                int idx3 = IDX_8D(k,j,g,r,i,j,a,r);
                                                int idx4 = IDX_8D(k,l,g,r,i,l,a,r);
                                                
                                                pi[idx1] += -im*U*G1y[idx2]*(G2xxxx[idx3]*kron_delta_bar(r,b) - kron_delta_bar(r,d)*G2xxxx[idx4]);
                                                pi[idx1] += -im*U*G1y[idx2]*(G2xyxy[idx3]*kron_delta_bar(r,b) - kron_delta_bar(r,d)*G2xyxy[idx4]);
                                            }
                                            int idx5 = IDX_8D(k,j,g,b,i,j,a,b);
                                            int idx6 = IDX_8D(k,l,g,d,i,l,a,d);
                                            pi[idx1] += -im*U*G1y[idx2]*(G2xyxy[idx5] - G2xyxy[idx6]);
                                        }
                                        else{
                                            if(l + n < Ns){
                                                for(int r = 0; r < Nb; ++r){
                                                    int idx8 = IDX_8D(k,l+n,g,r,i,l+n,a,r);
                                                    pi[idx1] += im*V[n-1]*G1y[idx2]*(G2xyxy[idx8] + G2xxxx[idx8]);
                                                }
                                            }
                                            if(l - n >= 0){
                                                for(int r = 0; r < Nb; ++r){
                                                    int idx9 = IDX_8D(k,l-n,g,r,i,l-n,a,r);
                                                    pi[idx1] += im*V[n-1]*G1y[idx2]*(G2xyxy[idx9] + G2xxxx[idx9]);
                                                }
                                            }
                                            
                                            if(j + n < Ns){
                                                for(int r = 0; r < Nb; ++r){
                                                    int idx10 = IDX_8D(k,j+n,g,r,i,j+n,a,r);
                                                    pi[idx1] += -im*V[n-1]*G1y[idx2]*(G2xyxy[idx10] + G2xxxx[idx10]);
                                                }
                                            }
                                            if(j - n >= 0){
                                                for(int r = 0; r < Nb; ++r){
                                                    int idx11 = IDX_8D(k,j-n,g,r,i,j-n,a,r);
                                                    pi[idx1] += -im*V[n-1]*G1y[idx2]*(G2xyxy[idx11] + G2xxxx[idx11]);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
#endif
