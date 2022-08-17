#include <mkl.h>
#include <thread>
#include <vector>
#include <complex>
using namespace std;

void hG_commutator(vector<complex<double>> &G,vector<complex<double>> &h, vector<complex<double>> &comm)
{
        //Function to calculate the commutator between the HF hamiltonian and the 1 particle greens function.
        char trans = 'N';
        complex<double> alpha,beta;
        alpha = 1.0;
        beta = 0.0;
        int N = int(sqrt(G.size()));
        vector<complex<double>> hG(N*N);
        vector<complex<double>> Gh(N*N);
        thread th11(cblas_zgemm,CblasRowMajor,CblasNoTrans,CblasNoTrans, N, N,N, &alpha,& *h.begin(),N,& *G.begin(),N,&beta,& *hG.begin(),N);
        thread th12(cblas_zgemm,CblasRowMajor,CblasNoTrans,CblasNoTrans, N, N,N, &alpha,& *G.begin(),N,& *h.begin(),N,&beta,& *Gh.begin(),N);
        th11.join();
        th12.join();
        for(int i = 0; i < N*N; ++i)
        {
                comm[i] = 0;
                comm[i] = hG[i] - Gh[i];
        }
}

void h2G2_commutator(vector<complex<double>> &G2, vector<complex<double>> &h1, vector<complex<double>> &comm, int Ns, int Nb)
{
        //Function that uses simplification given in corresponding note to calculate the commutator of the 2 particle Green's function with the tensor sum of two HF hamiltonians
        int Ns2 = int(pow(Ns,2));
	int Nb2 = int(pow(Nb,2));
	for(int i = 0; i < Ns; ++i)
        {
                for(int j = 0; j < Ns; ++j)
                {
                        for(int k = 0; k < Ns; ++k)
                        {
                                for(int l = 0; l < Ns;++l)
                                {
                                        for(int a = 0; a < Nb; ++a)
                                        {
                                                for(int b = 0; b< Nb; ++b)
                                                {
                                                        for(int g = 0; g<Nb; ++g)
                                                        {
                                                                for(int d = 0; d<Nb;++d)
                                                                {
                                                                        comm[(i+j*Ns)*Ns2*Nb2 + k+Ns*l + (a + b*Nb)*Ns2*Ns2*Nb2 + (g + d*Nb)*Ns2] = 0;

                                                                        for(int p = 0; p<Ns; ++p)
                                                                        {
                                                                                for(int r = 0; r<Nb; ++r)
                                                                                {
											comm[(i+j*Ns)*Ns2*Nb2 + k+Ns*l + (a + b*Nb)*Ns2*Ns2*Nb2 + (g + d*Nb)*Ns2] += h1[(j + b*Ns)*Nb*Ns + p + r*Ns]*G2[(i+p*Ns)*Ns2*Nb2 + k+ Ns*l + (a + r*Nb)*Ns2*Ns2*Nb2 + (g + d*Nb)*Ns2] ;
                                        						comm[(i+j*Ns)*Ns2*Nb2 + k+Ns*l + (a + b*Nb)*Ns2*Ns2*Nb2 + (g + d*Nb)*Ns2] += h1[(i + a*Ns)*Nb*Ns + p + r*Ns]*G2[(p+j*Ns)*Ns2*Nb2 + k+ Ns*l + (r + b*Nb)*Ns2*Ns2*Nb2 + (g + d*Nb)*Ns2] ;
                                        						comm[(i+j*Ns)*Ns2*Nb2 + k+Ns*l + (a + b*Nb)*Ns2*Ns2*Nb2 + (g + d*Nb)*Ns2] -= h1[(p + r*Ns)*Nb*Ns + l + d*Ns]*G2[(i+j*Ns)*Ns2*Nb2 + k+ Ns*p + (a + b*Nb)*Ns2*Ns2*Nb2 + (g + r*Nb)*Ns2] ;
                                        						comm[(i+j*Ns)*Ns2*Nb2 + k+Ns*l + (a + b*Nb)*Ns2*Ns2*Nb2 + (g + d*Nb)*Ns2] -= h1[(p + r*Ns)*Nb*Ns + k + g*Ns]*G2[(i+j*Ns)*Ns2*Nb2 + p+ Ns*l + (a + b*Nb)*Ns2*Ns2*Nb2 + (r + d*Nb)*Ns2] ;
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
