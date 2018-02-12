#pragma once

#include <fpdim3.h>

// geometrical and interacton parameters for spherocyls
struct SpherocylParameter
{
	double r;	// radius
	double t;	// interaction angle theta
	double d;	// interaction distance
};

class Spherocyl
{
	public:
		void setP(const fp::double3&);
		void setPN(const fp::double3&, const fp::double3&);
		void setN(const fp::double3&);
		void setS(const fp::double3&, const fp::double3&);

		fp::double3 getP();
		fp::double3 getN();
		fp::double3 getS1();
		fp::double3 getS2();

		unsigned int getID() const;
		void setID(unsigned int);

		void fixRoundingErrors();
		

	private:
		unsigned int ID;

		fp::double3 p;			// position, can be outside PBC
		fp::double3 n;			// normal direction, needs to be normalized
		fp::double3 s1;			// cylinder end position 1, can be outside PBC
		fp::double3 s2;			// cylinder end position 2, can be outside PBC

		void validateS();
		void validatePN();
};
