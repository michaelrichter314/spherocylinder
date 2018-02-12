#pragma once
#include <vector>
#include "mcmover.h"

#ifdef grid
#include "fpgrid.h"
#else
#include "kdtreespherocyl.h"
#endif

#include "spherocyl.h"
#include "cosinator.h"
#include "systemstate.h"
#include <fprandom.h>

using namespace std;
using namespace fp;

typedef vector<int> StaticCluster;
typedef SystemState<GridSize, Spherocyl, SpherocylParameter> SystemStateSpherocyl;

class Simulation
{
	public:
		Simulation(Randomizer*);

		void setState(const SystemState<GridSize, Spherocyl, SpherocylParameter>&);
		SystemState<GridSize, Spherocyl, SpherocylParameter> getState();

		double getE() const;
		double getV() const;
		void setBeta(const double&);

#ifdef NPT
		void setPressure(const double&);
#endif

		// debug stuff
		void checkSpace();

	protected:
#ifdef grid
		Grid space;							// grid holding the simulation box
#else
		KdTreeSpherocyl space;						// kd-tree holding the simulation box
#endif

		Randomizer* random = NULL;		// randomizer, externally provided by McMover

		vector<Spherocyl> spherocyl;	// array of all molecules in system
		unsigned int scN; 				// number of molecules in the system

		// spherocyl parameters
		SpherocylParameter spherocylParameters;	// initial spherocyl parameters, all other values are derived from this
		double scTipToTipDistSquared;
		double sc2RSquared;
		double scThetaToVFactor; //former scTfactor
		double scMaxInteractingCenterDistSquared;
		double scMaxInteractingEndDistSquared;
		double scMinGridSize;
		bool scMaxTwoInters;
		bool scAngleLimited;
		double scCosTheta;
#ifndef grid
		double intersectionBoxBorderSize;
		double interactionBoxBorderSize;
		double scMinSpaceDimension;
#endif

		void calcSpherocylParameters(const SpherocylParameter&);


		// thermodynamical parameters
		double E;						// system energy
		double beta;					// system inverse temperature
#ifdef NPT
		double P;						// system pressure
#endif


		// basic spherocylinder moves. move on grid, recalculate energy
		void move(Spherocyl&, const double3&);								// change position, update energy
		void rotate(Spherocyl&, const double3&, double, bool);				// rotate about one spherocyl end
		void rotate(Spherocyl&, const double3&, double, const double3&);	// rotate about specific point


		// intersection things
		pair<double3, double3> getIntersectionBox(Spherocyl&);
		bool intersect(Spherocyl&, Spherocyl&) const;
		bool intersects(Spherocyl&);


		// energy calculation things
		Cosinator cosinator;
		double calcInterTerm(const double&) const;		// calculates interaction term for given angle
#ifdef grid
		double calcInter(Spherocyl& a, Spherocyl& b) const;
#else
		pair<double3, double3> getInteractionBox(int);
		double calcInter(int, int);
#endif
		void calcE();									// sets E = energy. should only bew needed at beginning of simulation
		double calcEnergyContribution(Spherocyl& s);	// main energy calculation method


		// helpers for MC moves
		bool decide(const double&);						// decide to accept/reject a move according to energy difference
		// MC moves
		bool MCSingleGlobalTranslationRotation(const void*);
		bool MCSingleGlobalTranslation(const void*);
		bool MCSingleLocalTranslation(const void*);
		bool MCSingleRotation(const void*);
		bool MCSingleLocalRotation(const void*);
		bool MCSingleNormalTranslation(const void*);			// moves one spherocylinder in direction of its symmetry axis

		bool MCSinglePN(Spherocyl&, const double3&, const double3&); // helper function for moves of type that changes N and P for one spherocylinder
		
		bool MCClusterMove(const void*);

#ifdef NPT
		bool MCVolumeIso(const void*);
		bool MCVolumeAniso(const void*);
#endif

		// for cluster moves
		bool interact(Spherocyl&, Spherocyl&);
		bool interactsWithSomeone(Spherocyl&);
		StaticCluster getCluster(Spherocyl&);

		//friend class McMover;
		friend int main(int, char**);

#ifndef grid
		// kd tree traversal functions
		int globalId1;
		bool intersectsTraverse(int);
		bool interactsWithSomeoneTraverse(int);



#endif
};	
