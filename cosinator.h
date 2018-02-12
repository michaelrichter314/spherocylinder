#pragma once
#include <array>
#include <cmath>

class Cosinator
{
	public:
		Cosinator();
		double cos(const double&) const;
		double acos(const double&) const;

	private:
		static const unsigned int N = 65536; // 2^16
		const double PI = 3.14159265359;
		const double twoPi = 2.0*PI;

		std::array<double, N> cosValue;
		std::array<double, N> acosValue;
};
