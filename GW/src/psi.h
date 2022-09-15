#include "kron_delta.h"
#include <vector>
#include <complex>
using namespace std;


void Psi(vector<complex<double>> &psi, vector<complex<double>> &G1, int Ns, int Nb)
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
                                    int idx4 = IDX_4D(l,d,i,a);
                                    int idx5 = IDX_4D(k,g,j,b);
                                    psi[idx1] = (G1[idx2]*G1[idx3] - G1[idx4]*G1[idx5])*(kron_delta(i,j)*kron_delta_bar(a,b) - kron_delta(l,k)*kron_delta_bar(g,d));
                                    for(int r = 0; r<Nb; ++r){
                                        int idx6 = IDX_4D(j,r,i,a);
                                        int idx7 = IDX_4D(k,g,j,r);
                                        int idx8 = IDX_4D(l,r,i,a);
                                        int idx9 = IDX_4D(k,g,l,r);
                                        
                                        int idx10 = IDX_4D(i,r,j,b);
                                        int idx11 = IDX_4D(l,d,i,r);
                                        int idx12 = IDX_4D(k,r,j,b);
                                        int idx13 = IDX_4D(l,d,k,r);
                                        
                                        int idx14 = IDX_4D(j,r,i,a);
                                        int idx15 = IDX_4D(l,d,j,r);
                                        int idx16 = IDX_4D(k,r,i,a);
                                        int idx17 = IDX_4D(l,d,j,r);
                                        
                                        int idx18 = IDX_4D(i,r,j,b);
                                        int idx19 = IDX_4D(k,g,i,r);
                                        int idx20 = IDX_4D(l,r,j,b);
                                        int idx21 = IDX_4D(k,g,l,r);
                                        
                                        psi[idx1] += im*(G1[idx3]*(kron_delta_bar(r,b)*G1[idx6]*G1[idx7] - kron_delta_bar(r,d)*G1[idx8]*G1[idx9]));
                                        psi[idx1] += im*(G1[idx2]*(kron_delta_bar(r,a)*G1[idx10]*G1[idx11] - kron_delta_bar(r,g)*G1[idx12]*G1[idx13]));
                                        psi[idx1] -= im*(G1[idx5]*(kron_delta_bar(r,b)*G1[idx14]*G1[idx15] - kron_delta_bar(r,g)*G1[idx16]*G1[idx17]));
                                        psi[idx1] -= im*(G1[idx4]*(kron_delta_bar(r,a)*G1[idx18]*G1[idx19] - kron_delta_bar(r,d)*G1[idx20]*G1[idx21]));
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

