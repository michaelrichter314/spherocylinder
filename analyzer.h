#pragma once

#include "vtfloader.h"
#include <fpStatistics.h>

using namespace fp;
using namespace std;

class Analyzer
{
	public:
		Analyzer(string);
		void setTheta(double);
		bool goToNextTimestep();

		double getOrderParameter();
		double getAvNN();
		double getSmect();
		double getSmecticParameter();
		double getAssociationSpherocyl();
		double getAssociationEnd();

	private:
		double theta = 30;
		unique_ptr<VTFloader> loader;
		int timesteps;
		int timestep;

		SystemStateSpherocyl state;
		KdTreeSpherocyl space;

		void setState();

		double getNN(Spherocyl&);
		int countSmect(Spherocyl&);
		double3 getTrirector(Spherocyl&);

		// association fraction stuff
		pair<double3, double3> getInteractionBox(int);
		bool interact(int, int);
};
