#ifndef KRON_DELTA
#define KRON_DELTA

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

complex<double> kron_delta_lr(int x, int y, int n)
{
        //Kronecker delta complement for x and y
        return complex<double>(abs(x-y) == n);
}

#endif
