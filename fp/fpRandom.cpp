#include <fpRandom.h>

namespace fp
{
	bool randomSeeded = false;
#if fprandom == 1
	gsl_rng* randomNumberGenerator = NULL;
#endif
	// generate seed
	void randomSeed()
	{
            struct timeval time;
            gettimeofday(&time,NULL);
        #if fprandom == 1
            if (randomNumberGenerator != NULL)
            {
                // clear first
                gsl_rng_free (randomNumberGenerator);
            }
            
            //init
            randomNumberGenerator = gsl_rng_alloc (gsl_rng_mt19937);
            
            // generates a seed that changes every microsecond

            unsigned long int s = time.tv_sec*1000000 + time.tv_usec;
            
            gsl_rng_set (randomNumberGenerator, s);
            
            randomSeeded = true;
		#elif fprandom == 0
			// generates a seed that changes every millisecond

			srand((time.tv_sec * 1000) + (time.tv_usec / 1000)); 
            
            randomSeeded = true;
            
		#else
			#error Invalid fprandom
		#endif
	}
	
	// generates a double between 0 and 1, does not check if seeded
	inline double randomQuick()
	{
		#if fprandom == 0
			return (double)(rand()%RAND_MAX)/(double)RAND_MAX;
        #elif fprandom == 1
            return gsl_rng_uniform(randomNumberGenerator);
		#else
			#error Invalid fprandom
		#endif
	}
	
	// generates a double between 0 and 1
	inline double random()
	{
		if (!randomSeeded) {randomSeed();}
		
		return randomQuick();
	}
	
	// most efficient randomizer, works well for small numbers
	// does not check if random
	int randomQuick(int max)
	{
		#if fprandom == 0
			return rand()%max;
        #elif fprandom == 1
            return gsl_rng_get(randomNumberGenerator)%max;
		#else
			#error Invalid fprandom
		#endif
	}
	
	// Standard randomizer: checks if initialized, generates somewhat equally spaced random numbers
	int random(int max)
	{
		return ((int)(max*random()))%max;
	}
	
	// When things have to be guaranteed to be equally spaced, works well for large numbers
	int randomSlow(int max)
	{
		int toReturn = max;
		while(toReturn >= max)
		{
			#if fprandom == 0
				toReturn = rand();
#elif fprandom == 1
                toReturn = gsl_rng_get(randomNumberGenerator);
			#else
				#error Invalid fprandom
			#endif
		}
		
		return toReturn;
	}

    
    void randomQuit()
    {
#if fprandom == 1
        if (randomNumberGenerator != NULL)
        {
            gsl_rng_free (randomNumberGenerator);
        }
#endif
    }
}

