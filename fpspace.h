#pragma once

#include <fpdim3.h>
#include <iostream>
#include <cmath>
#include "spherocyl.h"

namespace fp
{
	// abstract class for space operations
	class Space
	{
		public:
			virtual void changeVolume(double3) = 0;
			virtual void insert(Spherocyl&) = 0;
			virtual void remove(Spherocyl&) = 0;

			double dist2PBC(const double3&, const double3&) const;
			double3 getPBC(const double3&) const;
			// displacement between 2 points that are inside PBC
			double3 dispPBC(const double3&, const double3&) const;
			// distance between two spherocyls of length 2 that are inside PBC
			double spherocylDist2PBC(const double3&, const double3&, const double3&, const double3&) const;

			const double3& getSpaceSize() const;

		protected:
			double3 spaceSize;
	};
}
