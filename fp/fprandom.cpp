#include <fprandom.h>

namespace fp
{
	Randomizer::Randomizer(int seed)
	{
		// creathe sphere buffer
		sphereBuffer = new double[3];

		struct timeval time;
		gettimeofday(&time,NULL);
		//init
		randomNumberGenerator = gsl_rng_alloc (gsl_rng_mt19937);

		// generates a seed that changes every microsecond
		unsigned long int s = time.tv_sec*1000000 + time.tv_usec + (unsigned long int)seed;

		gsl_rng_set (randomNumberGenerator, s);
	}

	Randomizer::~Randomizer()
	{
		gsl_rng_free (randomNumberGenerator);

		delete[] sphereBuffer;
	}

	double Randomizer::random()
	{
		return gsl_rng_uniform(randomNumberGenerator);
	}

	int Randomizer::random(int maxPlusOne)
	{
		return gsl_rng_get(randomNumberGenerator)%maxPlusOne;
	}

	double3 Randomizer::randomSpherical()
	{
		gsl_ran_dir_nd (randomNumberGenerator, 3, sphereBuffer);
		return double3{sphereBuffer[0], sphereBuffer[1], sphereBuffer[2]};
	}

	double3 Randomizer::randomSemispherical(double3 direction)
	{
		double3 temp;

		do
		{
			gsl_ran_dir_nd (randomNumberGenerator, 3, sphereBuffer);
			temp = double3{sphereBuffer[0], sphereBuffer[1], sphereBuffer[2]};
		} while (temp * direction < 0.0);

		return temp;
	}

	double3 Randomizer::randomSpherecap(double3 dir, double zmin, double zmax)
	{
		double z = zmin + (zmax-zmin)*random();
		return getWithFixedCos(dir, z);
	}

	double3 Randomizer::getWithFixedCos(double3 dir, double z)
	{
		double theta = random()*2.0*pi;
		double sq = sqrt(1.0-z*z);
		double x = sq*cos(theta);
		double y = sq*sin(theta);

		double3 toReturn{x,y,z};
		
		double3 n0 = {0.0,0.0,1.0};
		double3 n1 = normalize(dir);

		double cc=n1.z;
		if (fabs(1.0-cc) < 1e-6) {return toReturn;}
		if (fabs(-1.0-cc) < 1e-6) {return -toReturn;}

		double3 axis = cross(n0, n1);
		double angle = angleBetween(n1,n0);
		toReturn = rotate(toReturn, axis, angle);

		return toReturn;

	}

	int Randomizer::getSlot(double p, const std::vector<double>& v) const
	{
		int left = 0;
		int right = v.size()-1;

		if (p < v[left] || p > v[right])
		{
			std::cerr << p << " outside range " << v[left] << " " << v[right] << std::endl;
			return 0;
		}

		int i = (right+left)/2;
		//std::cout << p << std::endl;

		while (!(v[i+1] >= p && v[i] < p ))
		{
			//std::cout << p << " " << i << " " << left << " " << right << std::endl;
			i = (left + right)/2;
			if (v[i] < p) 
			{
				left = i;
			}
			else
			{
				right = i;
			}
		}
		//std::cout << i << " " << left << " " << right << std::endl;
		//std::cout << i << std::endl;
		return i;
	}

	int Randomizer::randomSlot(const std::vector<double>& v) //vec is integrated
	{
		double rnd = v[0] + random()*(v[v.size()-1] - v[0]);

		return getSlot(rnd, v);
	}

	double Randomizer::random(const std::vector<double>& v, int p)
	{
		double rMin = v[p];
		double rMax = v[p+1];
		/*if (p+1>=v.size())
		{
			std::cerr << "outside range" << std::endl;
		}*/
		return rMin + random()*(rMax - rMin);
	}

}
