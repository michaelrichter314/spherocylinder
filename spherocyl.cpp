#include "spherocyl.h"

void Spherocyl::setP(const fp::double3& newP)
{
	p = newP;
	validateS();
}

void Spherocyl::setPN(const fp::double3& newP, const fp::double3& newN)
{
	p = newP;
	n = newN;
	validateS();
}

void Spherocyl::setN(const fp::double3& newN)
{
	n = newN;
	validateS();
}

void Spherocyl::setS(const fp::double3& newS1, const fp::double3& newS2)
{
	s1 = newS1;
	s2 = newS2;

	validatePN();
}



fp::double3 Spherocyl::getP()
{
	return p;
}

fp::double3 Spherocyl::getN()
{
	return n;
}

fp::double3 Spherocyl::getS1()
{
	return s1;
}

fp::double3 Spherocyl::getS2()
{
	return s2;
}



unsigned int Spherocyl::getID() const
{
	return ID;
}

void Spherocyl::setID(unsigned int newID)
{
	ID = newID;
}



void Spherocyl::fixRoundingErrors()
{
	n = fp::normalize(n);
}



void Spherocyl::validateS()
{
	s1 = p + n;
	s2 = p - n;
}

void Spherocyl::validatePN()
{
	p = 0.5 * (s1 + s2);
	n = 0.5 * (s1 - s2);
}
