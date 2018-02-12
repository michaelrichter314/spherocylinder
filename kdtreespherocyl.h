#pragma once

#include <fpdim3.h>
#include <iostream>
#include <vector>
#include "fpspace.h"
#include "kdtree.h"

using namespace fp;
using namespace std;

class KdTreeSpherocyl : public Space
{
	public:
		void initialize(const double3&, Simulation*);
		void clear();

		// traverse whole tree until f is false
		void traverse(TraverseFunction);
		bool traverse(TraverseFunction, const double3&, const double3&);

		// return objects in certain region
		vector<int> getInBoxCOM(const double3&, const double3&);
		vector<int> getInBoxTip(const double3&, const double3&);


		// override stuff
		void insert(Spherocyl&);
		void remove(Spherocyl&);
		void changeVolume(double3);

	private:
		KdTree com;
		KdTree tip;
};

