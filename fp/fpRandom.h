#ifndef fprandom_h
#define fprandom_h

#ifndef fprandom
	#define fprandom 1
	// 0: std
	// 1: gsl
#endif

#if fprandom == 0
	#include <cstdlib>
	#include <time.h>
	#include <sys/time.h> 
#endif

#if fprandom == 1
	#include <time.h>
	#include <sys/time.h> 
	#include <gsl/gsl_rng.h>
#endif

namespace fp
{
	// generate seed
	void randomSeed();
	
	// generates a double between 0 and 1, does not check if seeded
	inline double randomQuick();
	
	// generates a double between 0 and 1
	double random();
	
	// most efficient randomizer, works well for small numbers
	// does not check if random
	int randomQuick(int max);
	
	// Standard randomizer: checks if initialized, generates somewhat equally spaced random numbers
	int random(int max);
	
	// When things have to be guaranteed to be equally spaced, works well for large numbers
	int randomSlow(int max);

    
    void randomQuit();
}

#endif
