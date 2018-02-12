#pragma once

#include <fpdim3.h>
#include <iostream>
#include <vector>
#include "fpspace.h"

using namespace fp;
using namespace std;

typedef double3 GridSize;
class Simulation;
using TraverseFunction = bool (Simulation::*)(int);

struct KdPoint
{
	double3 p;
	int N;
	int level;

	KdPoint* left = NULL;
	KdPoint* right = NULL;
	KdPoint* parent = NULL;
};

class KdTree
{
	public:
		void initialize(const double3&, Simulation*);
		~KdTree();

		void clear();

		// traverse whole tree until f is false
		void traverse(TraverseFunction);
		bool traverse(TraverseFunction, const double3&, const double3&);
		vector<int> getInBox(const double3&, const double3&);


		void insert(int, const double3&);
		void remove(int, const double3&);
		void changeVolume(double3);

	private:
		double3 spaceSize;
		KdPoint* root = NULL;

		KdPoint* find(int, const double3&);
		void remove(KdPoint*);
		KdPoint* findMin(KdPoint*, int);
		inline double getValue(KdPoint*, int);
		inline double getValue(const double3&, int);

		// temporary solution for a bug with S1, S2
		KdPoint* brute_force_find(int);
		KdPoint* brute_force_find(int, KdPoint*);

		// does not update parent pointers! only use to delete everything
		void removeSubtree(KdPoint*);

		bool traverse(KdPoint*, TraverseFunction);
		bool traverse(KdPoint*, TraverseFunction, const double3&, const double3&);
		void getInBox(KdPoint*, const double3&, const double3&, vector<int>&);

		// needed for traverse function
		Simulation* sim = NULL;
};

