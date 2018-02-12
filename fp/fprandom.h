#pragma once

#include <time.h>
#include <sys/time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>
#include <iostream>
#include <fpdim3.h>
#include <fpMath.h>

namespace fp
{
	class Randomizer
	{
		public:
			Randomizer(int = 0);
			~Randomizer();

			double random();
			int random(int);
			double3 randomSpherical();
			double3 randomSemispherical(double3);
			double3 randomSpherecap(double3, double, double);
			double3 getWithFixedCos(double3, double);
			int randomSlot(const std::vector<double>&); // picks a random number from a given distribution
			double random(const std::vector<double>&, int); // picks a random number in slot i
			int getSlot(double, const std::vector<double>&) const;
			

		private:
			gsl_rng* randomNumberGenerator;
			double* sphereBuffer;
	};
}
