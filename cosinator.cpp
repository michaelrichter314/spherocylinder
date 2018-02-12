#include "cosinator.h"

Cosinator::Cosinator()
{
	// fill arrays with life
	for (unsigned int i = 0; i < N; i++)
	{
		double x = (double)i/N*twoPi + 0.5/N*twoPi;
		cosValue[i] = std::cos(x);
	}

	for (unsigned int i = 0; i < N; i++)
	{
		double x = -1.0 + 2.0*(double)i/N + 1.0/N;
		acosValue[i] = std::acos(x);
	}
}

double Cosinator::cos(const double& x) const
{
	double xx = x - twoPi*std::floor(x/twoPi);

	return cosValue[(int)(xx/twoPi*N)];
}

double Cosinator::acos(const double& x) const
{
	if (x <= -1) {return acosValue[0];}
	if (x >= 1) {return acosValue[N-1];}

	return acosValue[(int)(N*(x+1.0)/2.0)];
}
	
