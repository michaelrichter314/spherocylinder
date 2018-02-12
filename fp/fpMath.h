#ifndef fpmath_h
#define fpmath_h

#ifndef fpnorm
	#define fpnorm euclid
#endif

#if fpnorm == default
	#define fpnorm euclid
#endif

#include <cmath>

namespace fp
{
	const double pi = 4*std::atan(1);
	
	
}

#endif
