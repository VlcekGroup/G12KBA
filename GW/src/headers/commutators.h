#include <mkl.h>
//#include <cblas.h>
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
#if		NOTHREADS
        cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, N, N,N, &alpha,& *h.begin(),N,& *G.begin(),N,&beta,& *hG.begin(),N);
        cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, N, N,N, &alpha,& *G.begin(),N,& *h.begin(),N,&beta,& *Gh.begin(),N);
#else
        thread th11(cblas_zgemm,CblasRowMajor,CblasNoTrans,CblasNoTrans, N, N,N, &alpha,& *h.begin(),N,& *G.begin(),N,&beta,& *hG.begin(),N);
        thread th12(cblas_zgemm,CblasRowMajor,CblasNoTrans,CblasNoTrans, N, N,N, &alpha,& *G.begin(),N,& *h.begin(),N,&beta,& *Gh.begin(),N);
        th11.join();
        th12.join();
#endif
        for(int i = 0; i < N*N; ++i)
        {
                comm[i] = hG[i] - Gh[i];
        }
}

void h2xyG2xyxy_commutator(vector<complex<double>> &G2xyxy,vector<complex<double>> &h1x, vector<complex<double>> &h1y, vector<complex<double>> &comm, int Ns, int Nb)
{
  //Function that uses simplification given in corresponding note to calculate the commutator of the 2 particle Green's function with the tensor sum of two HF hamiltonians
    int Ns2 = int(pow(Ns,2));
    int Nb2 = int(pow(Nb,2));

#define IDX_4D(x,y,z,t)  ((x) + (y)*Ns + ((z) + (t)*Ns)*Nb*Ns)
#define IDX_8D(x,y,z,t,m,n,q,w) ((x) + (y)*Ns + ((z) + (t)*Nb)*Ns2 + ((m)+(n)*Ns)*Ns2*Nb2 + ((q) + (w)*Nb)*Ns2*Ns2*Nb2)
#if NOTHREADS
  #pragma omp parallel for collapse(8)
#endif
    for(int b = 0; b < Nb; ++b){
        for(int a = 0; a < Nb; ++a){
            for(int j = 0; j < Ns; ++j){
                for(int i = 0; i < Ns; ++i){
                    for(int d = 0; d < Nb; ++d){
                        for(int g = 0; g < Nb; ++g){
                            for(int l = 0; l < Ns; ++l){
                                for(int k = 0; k < Ns; ++k){
                                    int idx1 = IDX_8D(k,l,g,d,i,j,a,b);
                                    comm[idx1] = 0.0;
                                    for(int r = 0; r<Nb; ++r){
                                        for(int p = 0; p<Ns; ++p){
                                            int idx2 = IDX_4D(p,r,j,b);
                                            int idx3 = IDX_4D(p,r,i,a);
                                            int idx4 = IDX_4D(l,d,p,r);
                                            int idx5 = IDX_4D(k,g,p,r);
                                            int idx6 = IDX_8D(k,l,g,d,i,p,a,r);
                                            int idx7 = IDX_8D(k,l,g,d,p,j,r,b);
                                            int idx8 = IDX_8D(k,p,g,r,i,j,a,b);
                                            int idx9 = IDX_8D(p,l,r,d,i,j,a,b);

                                            comm[idx1] += h1x[idx2]*G2xyxy[idx6];
                                            comm[idx1] += h1x[idx3]*G2xyxy[idx7];
                                            comm[idx1] -= h1y[idx4]*G2xyxy[idx8];
                                            comm[idx1] -= h1y[idx5]*G2xyxy[idx9];
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

