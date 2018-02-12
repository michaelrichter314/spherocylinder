#pragma once

#ifdef grid
#include "fpgrid.h"
#else
#include "kdtreespherocyl.h"
#endif
#include "spherocyl.h"
#include "systemstate.h"
#include <fprandom.h>
#include <iostream>
#include <vector>

using namespace std;
using namespace fp;

enum ICtype {ICrandom, ICsquarelattice, ICautomatic};

class ICmaker
{
	public:
		ICmaker(Randomizer*);
		ICmaker();
		~ICmaker();


		void setGridSize(const int3&);
		void setSpaceSize(const double3&);

		void setSpherocylParameters(const SpherocylParameter&);
		void setN(unsigned int);
		void setBeta(const double&);
#ifdef NPT
		void setP(const double&);
#endif

		void setType(const ICtype&);

		SystemState<GridSize, Spherocyl, SpherocylParameter> getState();


	private:
		ICtype type = ICrandom;
		int3 gridSize = int3{10,10,10};
		double3 spaceSize = double3{10.0, 10.0, 10.0};

		SpherocylParameter sp = SpherocylParameter {0.2, 30.0, 0.08};

		unsigned int N = 100;
		double beta;
#ifdef NPT
		double P;
#endif

		void make();
		SystemState<GridSize, Spherocyl, SpherocylParameter> state;
		vector<Spherocyl> spherocyl;

		bool hasBeenMade = false;

#ifdef grid
		Grid space;
#else
		KdTreeSpherocyl space;
#endif
		Randomizer* r = NULL;

		bool itsMyRandomizer = false;

		bool intersectsSomeone(Spherocyl&);
		bool intersect(Spherocyl&, Spherocyl&) const;
		pair<double3, double3> getIntersectionBox(Spherocyl&);
};

