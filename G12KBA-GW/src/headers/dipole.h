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

#ifndef dipole_h
#define dipole_h
#include <iostream>
#include <vector>
#include <complex>
using namespace std;

void dipole(vector<double> &p,vector<vector<complex<double>>> &G)
{
    int N = G[0].size();
	for(int i = 0; i < int(sqrt(N)/2); ++i){
        for(int j = 0; j < G.size(); ++j){
            		p[j]+=2*(sqrt(N)/2 - i - .5)*(G[j][i*(int(sqrt(N))+1)].imag() - G[j][N - i*(int(sqrt(N))+1)-1].imag())/sqrt(N);
		}
	}
}
#endif
