#include <vector>
#include <complex>
#include <iostream>
using namespace std;

void test_print(vector<complex<double>> &test_mat, int Ns, int Nb)
{
        for(int i = 0; i <Ns; ++i)
        {
                for(int j = 0; j < Ns; ++j)
                {
                        for(int a = 0; a < Nb; ++ a)
                        {
                                for(int b = 0;  b< Nb; ++b)
                                {
                                        cout<<"["<<i<<","<<j<<","<<a<<","<<b<<"]:"<<test_mat[(i + a*Ns)*Nb*Ns + (j+b*Ns)]<<"\n";
                                }
                        }
                }
        }
        cout<<"\n";
}
void test_print2(vector<complex<double>> &test_mat, int Ns,int Nb)
{
	int Ns2 = int(pow(Ns,2));
	int Nb2 = int(pow(Nb,2));
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
                                                        for(int g = 0; g<Nb; ++g)                                                                        {
                                                                for(int d = 0; d<Nb;++d)
                                                                {
                                        cout<<"["<<i<<","<<j<<","<<k<<","<<l<<":"<<a<<","<<b<<","<<g<<","<<d<<"]:"<< test_mat[(i+j*Ns)*Ns2*Nb2 + k+Ns*l + (a + b*Nb)*Ns2*Ns2*Nb2 + (g + d*Nb)*Ns2]<<"\n";
                                                                }
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }
}

