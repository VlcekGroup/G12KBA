#include <vector>
#include <complex>
using namespace std;


void collision_int(vector<complex<double>> &I, vector<complex<double>> &G2,int Ns, int Nb,double U)
{
        //Function to calculate the collision integral from the 2 particle greens function.  Implements equation 6 in corresponding note/equation 31 in Bonitz paper
        int Ns2 = int(pow(Ns,2));
	int Nb2 = int(pow(Nb,2));
	complex<double> im = {0,1};

	for(int i = 0; i<Ns; ++i)
        {
                for(int j = 0; j<Ns; ++j)
                {
                        for(int a = 0; a<Nb; ++a)
                        {
                                for(int b = 0; b<Nb; ++b)
                                {
                                        I[(i + a*Ns)*Ns*Nb + j + b*Ns] = 0;
                                        for(int m = 0; m<Nb;++m)
                                        {
                                                if(m != a)
                                                {
                                                        I[(i+a*Ns)*Nb*Ns + j + b*Ns] += -im*U*G2[(i + i*Ns)*Ns2*Nb2 + j + Ns*i + (a + m*Nb)*Ns2*Ns2*Nb2 + (b + m*Nb)*Ns2];
                                                }
                                        }
                                }
                        }
                }
        }
}
