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


#ifndef kron_delta_h
#define kron_delta_h

#include <complex>
using namespace std;
complex<double> kron_delta(int x, int y)
{
        //Kronecker delta complement for x and y
        return complex<double>(x==y);
}


complex<double> kron_delta_bar(int x, int y)
{
        //Kronecker delta complement for x and y
        return complex<double>(x!=y);
}

#endif
