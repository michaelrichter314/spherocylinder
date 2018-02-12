#include "fpspace.h"

namespace fp
{

double Space::dist2PBC(const double3& a, const double3& b) const
	{
		const auto& size = spaceSize;

		double3 up = {std::fmax(a.x,b.x), std::fmax(a.y,b.y), std::fmax(a.z,b.z)};
		double3 lo = {std::fmin(a.x,b.x), std::fmin(a.y,b.y), std::fmin(a.z,b.z)};

		double dx = std::fmin(up.x - lo.x, lo.x + size.x - up.x);
		double dy = std::fmin(up.y - lo.y, lo.y + size.y - up.y);
		double dz = std::fmin(up.z - lo.z, lo.z + size.z - up.z);

		return dx*dx + dy*dy + dz*dz;
	}

	double3 Space::getPBC(const double3& a) const
	{
		double dx = a.x - spaceSize.x*std::floor(a.x/spaceSize.x);
		double dy = a.y - spaceSize.y*std::floor(a.y/spaceSize.y);
		double dz = a.z - spaceSize.z*std::floor(a.z/spaceSize.z);

		return double3{dx,dy,dz};
	}

	double3 Space::dispPBC(const double3& aa, const double3& bb) const
	{
		double3 rij = bb - aa;

		double3 res;

		double FX;
		FX = std::fabs(rij.x);
		if (FX < spaceSize.x-FX)
		{
			res.x = rij.x;
		}
		else
		{
			res.x = rij.x - ((rij.x > 0)?spaceSize.x:-spaceSize.x);
		}

		FX = std::fabs(rij.y);
		if (FX < spaceSize.y-FX)
		{
			res.y = rij.y;
		}
		else
		{
			res.y = rij.y - ((rij.y > 0)?spaceSize.y:-spaceSize.y);
		}

		FX = std::fabs(rij.z);
		if (FX < spaceSize.z-FX)
		{
			res.z = rij.z;
		}
		else
		{
			res.z = rij.z - ((rij.z > 0)?spaceSize.z:-spaceSize.z);
		}

		return res;

	}

	double Space::spherocylDist2PBC(const double3& p1, const double3& w1, const double3& p2, const double3& w2) const
	{
		double3 min_rij = dispPBC(p1, p2);

		double xla, xmu;

		double rr = norm2(min_rij);
		double rw1 = min_rij*w1;
		double rw2 = min_rij*w2;
		double w1w2 = w1*w2;
		double cc = 1.0 - w1w2*w1w2;

		const double cutoff = 1e-6;

		if (cc < cutoff)
		{
			// if roughly parallel
			if (rw1 == 0.0 || rw2 == 0.0)
			{
				// if center of mass next to each other
				return rr;
			}
			else
			{
				xla = rw1/2;
				xmu = -rw2/2;
			}
		}
		else
		{
			// step 1
			xla = (rw1-w1w2*rw2) / cc;
			xmu = (-rw2+w1w2*rw1) / cc;
		}

		// step2
		if ( fabs(xla) > 1.0 || fabs(xmu) > 1.0 )
		{
			// step 3-7

			if( fabs(xla)-1.0 > fabs(xmu)-1.0 )
			{
				xla = (xla<0)?-1.0:1.0;
				xmu = xla*w1w2-rw2;
				if ( fabs(xmu)>1.0 )
				{
					xmu = (xmu<0)?-1.0:1.0;
				}
			}
			else
			{
				xmu = (xmu<0)?-1.0:1.0;;
				xla = xmu*w1w2+rw1;
				if ( fabs(xla)>1.0 )
				{
					xla = (xla<0)?-1.0:1.0;
				}
			}
		}

		// step 8
		return rr+xla*xla+xmu*xmu + 2*(xmu*rw2 -xla*(rw1+xmu*w1w2));
	}
	
	const double3& Space::getSpaceSize() const
	{
		return spaceSize;
	}
}
