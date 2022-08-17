#include <vector>
#include <complex>
using namespace std;
void G1_step(vector<complex<double>> &h1G1_comm, vector<complex<double>> &I,vector<complex<double>> &K,vector<complex<double>> &RK_result, double w,int Ns,int Nb)
{
        //Function takes a partial Runge Kutta time step for the full matrix of G1 and stores the result in K for use in the next RK step as well as in RK_result which is to be used
        //in updating G1 at the end of each full RK time step 
        complex<double> im = {0,1};
	for(int i = 0; i < Ns; ++i)
        {
                for(int j = 0; j< Ns; ++j)
                {
                        for(int a = 0; a<Nb;++a)
                        {
                                for(int b = 0; b<Nb;++b)
                                {
                                        K[(i + a*Ns)*Ns*Nb + j + b*Ns] = -im*(h1G1_comm[(i+Ns*a)*Nb*Ns + j + b*Ns] + (I[(i+Ns*a)*Nb*Ns + j + b*Ns]) + conj(I[(j+Ns*b)*Nb*Ns + i + a*Ns]));
                                        RK_result[(i + a*Ns)*Nb*Ns + j + b*Ns] += w*K[(i + a*Ns)*Nb*Ns + j + b*Ns];
                                }
                        }
                }
        }
}






void G2_step(vector<complex<double>> &h2G2_comm, vector<complex<double>> &psi,vector<complex<double>> &K,vector<complex<double>> &RK_result, double w,int Ns,int Nb,double U)
{
        //Function implements same time stepping as G1_step but for the 2 particle green's function
	int Ns2 = int(pow(Ns,2));
	int Nb2 = int(pow(Nb,2));
        complex<double> im = {0,1};
	for(int i = 0; i < Ns; ++i)
        {
                for(int j = 0; j< Ns; ++j)
                {
                        for(int k = 0; k < Ns; ++k)
                        {
                                for(int l = 0; l< Ns; ++l)
                                {
                                        for(int a = 0; a < Nb; ++a)
                                        {
                                                for(int b = 0; b< Nb; ++b)
                                                {
                                                        for(int g = 0; g < Nb; ++g)
                                                        {
                                                                for(int d = 0; d< Nb; ++d)
                                                                {
									K[(i + j*Ns)*Ns2*Nb2 + k + Ns*l + (a + b*Nb)*Ns2*Ns2*Nb2 + (g + d*Nb)*Ns2] = -im*(h2G2_comm[(i+j*Ns)*Ns2*Nb2 + k+Ns*l + (a + b*Nb)*Ns2*Ns2*Nb2 + (g + d*Nb)*Ns2] +  U*psi[(i+j*Ns)*Ns2*Nb2 + k+Ns*l + (a + b*Nb)*Ns2*Ns2*Nb2 + (g + d*Nb)*Ns2]);
                                                			RK_result[(i+j*Ns)*Ns2*Nb2 + k+Ns*l + (a + b*Nb)*Ns2*Ns2*Nb2 + (g + d*Nb)*Ns2] += w*K[(i+j*Ns)*Ns2*Nb2 + k+Ns*l + (a + b*Nb)*Ns2*Ns2*Nb2 + (g + d*Nb)*Ns2];
                                                                }
                                                        }
                                                }
                                        }
                                }
                        }
                 }
        }
}

