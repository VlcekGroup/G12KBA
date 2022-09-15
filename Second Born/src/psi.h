#include <vector>
#include <complex>
using namespace std;
complex<double> kron_delta(int x, int y)
{
        //Function that returns the Kronecker delta for x and y
        return complex<double>(x == y);
}

complex<double> kron_delta_bar(int x, int y)
{
        //Kronecker delta complement for x and y
        return complex<double>(x!=y);
}

void Psi(vector<complex<double>> &psi, vector<complex<double>> &G1, int Ns, int Nb)
{
	//Function that calculates what is called the 2 particle occupations in G1-G2 paper by Bonitz.  Implements equation 7 in corresponding note/unlabelled equation after equation 21 in Bonitz paper
	int Ns2 = int(pow(Ns,2));
	int Nb2 = int(pow(Nb,2));
	complex<double> im = {0,1};
#if NOTHREADS
  #pragma omp parallel for collapse(8)
#endif
	for(int i = 0; i < Ns; ++i)
	{
	        for(int j = 0; j< Ns; ++j)
	        {
	                for(int k = 0; k< Ns; ++k)
	                {
	                        for(int l = 0; l<Ns;++l)
	                        {
	                                for(int a = 0; a< Nb; ++a)
	                                {
	                                        for(int b = 0; b< Nb; ++b)
	                                        {
	                                                for(int g = 0; g<Nb; ++g)
	                                                {
	                                                        for(int d = 0; d<Nb;++d)
	                                                        {
                                                                psi[(i+j*Ns)*Ns2*Nb2 + k+Ns*l + (a + b*Nb)*Ns2*Ns2*Nb2 + (g + d*Nb)*Ns2] = (G1[(i + a*Ns)*Nb*Ns + (k + g*Ns)]*G1[(j+b*Ns)*Nb*Ns + (l + d*Ns)])*(kron_delta(i,j)*kron_delta_bar(a,b) - kron_delta(l,k)*kron_delta_bar(g,d));
	//                                                                                      psi[(i+j*Ns)*Ns2*Nb2 + k+Ns*l + (a + b*Nb)*Ns2*Ns2*Nb2 + (g + d*Nb)*Ns2] = (G1[(i + a*Ns)*Nb*Ns + (k + g*Ns)]*G1[(j+b*Ns)*Nb*Ns + (l + d*Ns)] - G1[(i$
	                                                                for(int r = 0; r<Nb; ++r)
	                                                                {
                                                                        psi[(i+j*Ns)*Ns2*Nb2 + k+Ns*l + (a + b*Nb)*Ns2*Ns2*Nb2 + (g + d*Nb)*Ns2] += im*(G1[(j + b*Ns)*Nb*Ns + (l + d*Ns)]*(kron_delta_bar(r,b)*G1[(i + a*Ns)*Nb*Ns + (j + r*Ns)]*G1[(j + r*Ns)*Nb*Ns + (k + g*Ns)] - (kron_delta_bar(r,d))*G1[(i + a*Ns)*Nb*Ns + (l + r*Ns)]*G1[(l + r*Ns)*Nb*Ns + (k + g*Ns)]));
                                                                        psi[(i+j*Ns)*Ns2*Nb2 + k+Ns*l + (a + b*Nb)*Ns2*Ns2*Nb2 + (g + d*Nb)*Ns2] += im*(G1[(i + a*Ns)*Nb*Ns + (k + g*Ns)]*((kron_delta_bar(r,a))*G1[(j + b*Ns)*Nb*Ns + (i + r*Ns)]*G1[(i + r*Ns)*Nb*Ns + (l + d*Ns)] - (kron_delta_bar(r,g))*G1[(j + b*Ns)*Nb*Ns + (k + r*Ns)]*G1[(k + r*Ns)*Nb*Ns + (l + d*Ns)]));

	                                                                        //psi[(i+j*Ns)*Ns2*Nb2 + k+Ns*l + (a + b*Nb)*Ns2*Ns2*Nb2 + (g + d*Nb)*Ns2] -= im*(G1[(i + a*Ns)*Nb*Ns + (l + d*Ns)]*(kron_delta_bar(r,a)*G1[($
	                                                                        //psi[(i+j*Ns)*Ns2*Nb2 + k+Ns*l + (a + b*Nb)*Ns2*Ns2*Nb2 + (g + d*Nb)*Ns2] -= im*(G1[(j + b*Ns)*Nb*Ns + (k + g*Ns)]*(kron_delta_bar(r,b)*G1[($
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

