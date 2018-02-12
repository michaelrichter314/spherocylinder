#include "fpdim3.h"

namespace fp
{
	

	double3 operator+ (const double3& lhs, const int3& rhs)
	{
		double3 toReturn = lhs;
		toReturn.x += (double)rhs.x;
		toReturn.y += (double)rhs.y;
		toReturn.z += (double)rhs.z;

		return toReturn;
	}

	double3 operator+ (const int3& lhs, const double3& rhs)
	{
		double3 toReturn = rhs;
		toReturn.x += (double)lhs.x;
		toReturn.y += (double)lhs.y;
		toReturn.z += (double)lhs.z;

		return toReturn;
	}

	double3 operator+= (double3& lhs, const int3& rhs)
	{
		lhs.x += (double)rhs.x;
		lhs.y += (double)rhs.y;
		lhs.z += (double)rhs.z;

		return lhs;
	}

	double3 operator*(const int3& lhs, const double& rhs)
	{
		double3 toReturn;
		toReturn.x = lhs.x*rhs;
		toReturn.y = lhs.y*rhs;
		toReturn.z = lhs.z*rhs;
		return toReturn;
	}

	double3 operator*(const double& lhs, const int3& rhs)
	{
		double3 toReturn;
		toReturn.x = rhs.x*lhs;
		toReturn.y = rhs.y*lhs;
		toReturn.z = rhs.z*lhs;
		return toReturn;
	}

	double3 normalize(const double3& in)
	{
		return in * (1.0 / norm(in));
	}

	double angleBetween(const double3& a, const double3& b)
	{
		return 180.0/pi * acos( normalize(a) * normalize(b) );
	}
	
	double radBetween(const double3& a, const double3& b)
	{
		return acos( normalize(a) * normalize(b) );
	}

	// Rodrigues' rotation formula
	double3 rotate(const double3& v, const double3& a, double t)
	{
		double3 k = normalize(a);
		double th = t * 3.14159265359 / 180.0;

		return v * cos(th) + cross(k, v)*sin(th) + k * (k*v) * (1.0 - cos(th));
	}

	int3 round(const double3& d)
	{
		return {(int)(std::round(d.x)), (int)(std::round(d.y)), (int)(std::round(d.z))};
	}




}

