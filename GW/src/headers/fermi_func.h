#ifndef FERMI_FUNC_H
#define FERMI_FUNC_H
#include <cmath>
double fermi(double t, double t0,double tau)
{
	return 1.0/(1+exp((t0-t)/tau));
}

double erf(double t,double t0,double beta)
{
	return (erf(beta*(t - t0))+1.0)/2.0;
}

double opt(double t,double N, double t0)
{
return (1.0 + tan((2.0*t/t0-1)*atan(sqrt(N-1)))/sqrt(N-1))/2.0;

}
#endif
