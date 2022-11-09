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

#ifndef psi_h
#define psi_h

#include "kron_delta.h"
#include <vector>
#include <complex>
using namespace std;

void Psi_xyxy(vector<complex<double>> &psi, vector<complex<double>> &G1, int Ns, int Nb,double U, vector<double> &V)
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
                                    int idx2 = IDX_4D(k,g,i,a);
                                    int idx3 = IDX_4D(l,d,j,b);
                                    psi[idx1] = U*G1[idx2]*G1[idx3]*(kron_delta(i,j)*kron_delta_bar(a,b) - kron_delta(l,k)*kron_delta_bar(g,d));
                                    for(int n = 0; n < V.size() + 1; ++n){
                                        if(n==0){
                                            for(int r = 0; r<Nb; ++r){
                                                int idx4 = IDX_4D(j,r,i,a);
                                                int idx5 = IDX_4D(k,g,j,r);
                                                int idx6 = IDX_4D(l,r,i,a);
                                                int idx7 = IDX_4D(k,g,l,r);
                                                int idx8 = IDX_4D(i,r,j,b);
                                                int idx9 = IDX_4D(l,d,i,r);
                                                int idx10 = IDX_4D(k,r,j,b);
                                                int idx11 = IDX_4D(l,d,k,r);
                                                
                                                psi[idx1] += im*U*(G1[idx2]*(kron_delta_bar(r,a)*G1[idx8]*G1[idx9] - kron_delta_bar(r,g)*G1[idx10]*G1[idx11]));
                                                psi[idx1] += im*U*(G1[idx3]*(kron_delta_bar(r,b)*G1[idx4]*G1[idx5] - kron_delta_bar(r,d)*G1[idx6]*G1[idx7]));
                                            }
                                            int idx12 = IDX_4D(l,d,i,a);
                                            int idx13 = IDX_4D(k,g,j,b);
                                            int idx14 = IDX_4D(l,d,k,g);
                                            int idx15 = IDX_4D(i,a,j,b);
                                            int idx16 = IDX_4D(k,g,l,d);
                                            int idx17 = IDX_4D(j,b,i,a);
                                            
                                            psi[idx1]+= -G1[idx2]*U*G1[idx3]*(kron_delta(k,l)*kron_delta(d,g) - kron_delta(i,j)*kron_delta(a,b));
                                            psi[idx1]+= -im*G1[idx2]*U*(G1[idx13]*G1[idx14] - G1[idx15]*G1[idx12]);
                                            psi[idx1]+= -im*G1[idx3]*U*(G1[idx12]*G1[idx16] - G1[idx17]*G1[idx13]);
                                        }
                                        else{
                                            
                                            if(abs(k - l) == n){
                                                psi[idx1] += -V[n-1]*G1[idx2]*G1[idx3];
                                            }
                                            if(abs(i - j) == n){
                                                psi[idx1] += V[n-1]*G1[idx2]*G1[idx3];
                                            }
                                            
                                    
                                            if(l + n < Ns){
                                                for(int r = 0; r < Nb; ++r){
                                                    int idx18 = IDX_4D(l+n,r,i,a);
                                                    int idx19 = IDX_4D(k,g,l+n,r);
                                                    psi[idx1] += -im*V[n-1]*G1[idx18]*G1[idx19]*G1[idx3];
                                                }
                                            }
                                            if(l - n >= 0){
                                                for(int r = 0; r < Nb; ++r){
                                                    int idx20 = IDX_4D(l-n,r,i,a);
                                                    int idx21 = IDX_4D(k,g,l-n,r);
                                                    psi[idx1] += -im*V[n-1]*G1[idx20]*G1[idx21]*G1[idx3];
                                                }
                                            }
                                            
                                            if(j + n < Ns){
                                                for(int r = 0; r < Nb; ++r){
                                                    int idx22 = IDX_4D(j+n,r,i,a);
                                                    int idx23 = IDX_4D(k,g,j+n,r);
                                                    psi[idx1] += im*V[n-1]*G1[idx22]*G1[idx23]*G1[idx3];
                                                }
                                            }
                                            if(j - n >= 0){
                                                for(int r = 0; r < Nb; ++r){
                                                    int idx24 = IDX_4D(j-n,r,i,a);
                                                    int idx25 = IDX_4D(k,g,j-n,r);
                                                    psi[idx1] += im*V[n-1]*G1[idx24]*G1[idx25]*G1[idx3];
                                                }
                                            }
                                            if(k + n < Ns){
                                                for(int r = 0; r < Nb; ++r){
                                                    int idx26 = IDX_4D(k+n,r,j,b);
                                                    int idx27 = IDX_4D(l,d,k+n,r);
                                                    psi[idx1] += -im*V[n-1]*G1[idx26]*G1[idx27]*G1[idx2];
                                                }
                                            }
                                            if(k - n >= 0){
                                                for(int r = 0; r < Nb; ++r){
                                                    int idx28 = IDX_4D(k-n,r,j,b);
                                                    int idx29 = IDX_4D(l,d,k-n,r);
                                                    psi[idx1] += -im*V[n-1]*G1[idx28]*G1[idx29]*G1[idx2];
                                                }
                                            }
                                            
                                            if(i + n < Ns){
                                                for(int r = 0; r < Nb; ++r){
                                                    int idx30 = IDX_4D(i+n,r,j,b);
                                                    int idx31 = IDX_4D(l,d,i+n,r);
                                                    psi[idx1] += im*V[n-1]*G1[idx30]*G1[idx31]*G1[idx2];
                                                }
                                            }
                                            if(i - n >= 0){
                                                for(int r = 0; r < Nb; ++r){
                                                    int idx32 = IDX_4D(i-n,r,j,b);
                                                    int idx33 = IDX_4D(l,d,i-n,r);
                                                    psi[idx1] += im*V[n-1]*G1[idx32]*G1[idx33]*G1[idx2];
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


void Psi_xxxx(vector<complex<double>> &psi, vector<complex<double>> &G1, int Ns, int Nb,double U, vector<double> &V)
{
    //Function that calculates what is called the 2 particle occupations in G1-G2 paper by Bonitz.  Implements equation 7 in corresponding note/unlabelled equation after equation 21 in Bonitz paper
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
                                    int idx2 = IDX_4D(k,g,i,a);
                                    int idx3 = IDX_4D(l,d,j,b);
                                    
                                    psi[idx1] = U*G1[idx2]*G1[idx3]*(kron_delta(i,j)*kron_delta_bar(a,b) - kron_delta(l,k)*kron_delta_bar(g,d));
                                    for(int n = 0; n < V.size() + 1; ++n){
                                        if(n==0){
                                            for(int r = 0; r<Nb; ++r){
                                                int idx4 = IDX_4D(j,r,i,a);
                                                int idx5 = IDX_4D(k,g,j,r);
                                                int idx6 = IDX_4D(l,r,i,a);
                                                int idx7 = IDX_4D(k,g,l,r);
                                                
                                                int idx8 = IDX_4D(i,r,j,b);
                                                int idx9 = IDX_4D(l,d,i,r);
                                                int idx10 = IDX_4D(k,r,j,b);
                                                int idx11 = IDX_4D(l,d,k,r);
                                                
                                                psi[idx1] += im*(G1[idx3]*U*(kron_delta_bar(r,b)*G1[idx4]*G1[idx5] - kron_delta_bar(r,d)*G1[idx6]*G1[idx7]));
                                                psi[idx1] += im*(G1[idx2]*U*(kron_delta_bar(r,a)*G1[idx8]*G1[idx9] - kron_delta_bar(r,g)*G1[idx10]*G1[idx11]));
                                            }
                                        }
                                        else{
                                            
                                            if(abs(k - l) == n){
                                                psi[idx1] += -V[n-1]*G1[idx2]*G1[idx3];
                                            }
                                            if(abs(i - j) == n){
                                                psi[idx1] += V[n-1]*G1[idx2]*G1[idx3];
                                            }
                                            
                                            
                                            if(l + n < Ns){
                                                for(int r = 0; r < Nb; ++r){
                                                    int idx18 = IDX_4D(l+n,r,i,a);
                                                    int idx19 = IDX_4D(k,g,l+n,r);
                                                    psi[idx1] += -im*V[n-1]*G1[idx18]*G1[idx19]*G1[idx3];
                                                }
                                            }
                                            if(l - n >= 0){
                                                for(int r = 0; r < Nb; ++r){
                                                    int idx20 = IDX_4D(l-n,r,i,a);
                                                    int idx21 = IDX_4D(k,g,l-n,r);
                                                    psi[idx1] += -im*V[n-1]*G1[idx20]*G1[idx21]*G1[idx3];
                                                }
                                            }
                                            
                                            if(j + n < Ns){
                                                for(int r = 0; r < Nb; ++r){
                                                    int idx22 = IDX_4D(j+n,r,i,a);
                                                    int idx23 = IDX_4D(k,g,j+n,r);
                                                    psi[idx1] += im*V[n-1]*G1[idx22]*G1[idx23]*G1[idx3];
                                                }
                                            }
                                            if(j - n >= 0){
                                                for(int r = 0; r < Nb; ++r){
                                                    int idx24 = IDX_4D(j-n,r,i,a);
                                                    int idx25 = IDX_4D(k,g,j-n,r);
                                                    psi[idx1] += im*V[n-1]*G1[idx24]*G1[idx25]*G1[idx3];
                                                }
                                            }
                                            if(k + n < Ns){
                                                for(int r = 0; r < Nb; ++r){
                                                    int idx26 = IDX_4D(k+n,r,j,b);
                                                    int idx27 = IDX_4D(l,d,k+n,r);
                                                    psi[idx1] += -im*V[n-1]*G1[idx26]*G1[idx27]*G1[idx2];
                                                }
                                            }
                                            if(k - n >= 0){
                                                for(int r = 0; r < Nb; ++r){
                                                    int idx28 = IDX_4D(k-n,r,j,b);
                                                    int idx29 = IDX_4D(l,d,k-n,r);
                                                    psi[idx1] += -im*V[n-1]*G1[idx28]*G1[idx29]*G1[idx2];
                                                }
                                            }
                                            
                                            if(i + n < Ns){
                                                for(int r = 0; r < Nb; ++r){
                                                    int idx30 = IDX_4D(i+n,r,j,b);
                                                    int idx31 = IDX_4D(l,d,i+n,r);
                                                    psi[idx1] += im*V[n-1]*G1[idx30]*G1[idx31]*G1[idx2];
                                                }
                                            }
                                            if(i - n >= 0){
                                                for(int r = 0; r < Nb; ++r){
                                                    int idx32 = IDX_4D(i-n,r,j,b);
                                                    int idx33 = IDX_4D(l,d,i-n,r);
                                                    psi[idx1] += im*V[n-1]*G1[idx32]*G1[idx33]*G1[idx2];
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
