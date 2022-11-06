#include <iostream>
#include <vector>
#include <complex>
using namespace std;
void dipole(vector<double> &p,vector<vector<complex<double>>> &G)
{
	for(int i = 0; i < int(sqrt(G[0].size())/2); ++i){
        	int N = G[0].size();
        	for(int j = 0; j < G.size(); ++j){
            		p[j]+=2*(sqrt(N)/2 - i - .5)*(G[j][i*(int(sqrt(G[0].size()))+1)].imag() - G[j][N - i*(int(sqrt(G[0].size()))+1)-1].imag())/sqrt(G[0].size());
		}
	}
}
