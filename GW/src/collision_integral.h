#include <vector>
#include <complex>
using namespace std;


void collision_int(vector<complex<double>> &I, vector<complex<double>> &G2,int Ns, int Nb,double U)
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
                    I[idx1] = 0;
                    for(int m = 0; m<Nb;++m){
                        if(m != a){
                            int idx2 = IDX_8D(j,i,b,m,i,i,a,m);
                            I[idx1] += -im*U*G2[idx2];
                        }
                    }
                }
            }
        }
    }
}


