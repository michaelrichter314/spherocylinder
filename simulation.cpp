#include "simulation.h"

Simulation::Simulation(Randomizer* r) : random(r)
{
}

void Simulation::calcSpherocylParameters(const SpherocylParameter& p)
{
	const double& r = p.r;
	const double& t = p.t;
	const double& d = p.d;
	// r=radius
	// t=theta in degrees
	// d=interaction distance
	//
	
	spherocylParameters = p;
	
	const double PI = 3.14159265359;
	
	scTipToTipDistSquared = (2*r+2)*(2*r+2);
	sc2RSquared = 4*r*r;
	scThetaToVFactor = 180.0/t;
	scMaxTwoInters = (4*r+2*d < 2.0); //TODO: use theta
	scMaxInteractingCenterDistSquared = (2+2*r+d)*(2+2*r+d);
	scMaxInteractingEndDistSquared = (2*r+d)*(2*r+d);
	scAngleLimited = ( t < 180.0 );
	scCosTheta = cos(t*PI/180.0);
	scMinGridSize = 2.0 + 2*r + d;
#ifndef grid
	intersectionBoxBorderSize = 1.0 + 2*r;
	interactionBoxBorderSize = 2*r+d;
	scMinSpaceDimension = fmax(4.0 + 4.0*r, 2.0*d+4.0*r);
#endif
	
}
void Simulation::move(Spherocyl& s, const double3& d)
{
	double3 oldPos = s.getP();
	double3 newPos = space.getPBC(oldPos + d);

	// update spherocylinder info
	space.remove(s);
	s.setP(newPos);

	// update grid
	space.insert(s);
}

void Simulation::rotate(Spherocyl& s, const double3& axis, double angle, bool orig)
{

	double3 oldN = s.getN();
	double3 newN = fp::rotate(oldN, axis, angle);

	double3 s1disp = newN - oldN;
	double3 oldPos = s.getP();
	double3 newPos;

	if (orig)
	{
		newPos = space.getPBC(oldPos - s1disp);
	}
	else
	{
		newPos = space.getPBC(oldPos + s1disp);
	}

	// update spherocylinder info
	space.remove(s);
	s.setPN(newPos, newN);

	// update grid
	space.insert(s);
}

void Simulation::rotate(Spherocyl& s, const double3& axis, double angle, const double3& orig)
{

	double3 oldN = s.getN();
	double3 newN = fp::rotate(oldN, axis, angle);

	double3 oldPos = s.getP();
	double3 newPos = space.getPBC(orig + fp::rotate(space.dispPBC(orig,oldPos), axis, angle));

	// update spherocylinder info
	space.remove(s);
	s.setPN(newPos, newN);

	// update grid
	space.insert(s);
}



bool Simulation::intersect(Spherocyl& a, Spherocyl& b) const
{
	double3 p1 = a.getP();
	double3 p2 = b.getP();

	// if they are close, calculate actual distance
	if (space.dist2PBC(p1, p2) <= scTipToTipDistSquared)
	{
		// check actual distance
		if (space.spherocylDist2PBC(p1, a.getN(), p2, b.getN()) <= sc2RSquared)
		{
			return true;
		}
	}
	return false;
}

bool Simulation::intersects(Spherocyl& s) // not const since grid iterator changes
{
#ifdef grid
	int id1 = s.getID();
	space.getIDsNear(s.getP());

	while (space.gotMoreIDs())
	{
		int id2 = space.getNextID();
		if (id1 != id2 && intersect(s, spherocyl[id2]))
		{
			return true;
		}
	}

	return false;
#else
	// Kd-tree
	auto box = getIntersectionBox(s);

	globalId1 = s.getID();
	return (!space.traverse(&Simulation::intersectsTraverse, box.first, box.second));
#endif
}

double Simulation::calcInterTerm(const double& c) const
{
	//double angle = cosinator.acos(c);
	//double cc = cosinator.cos(angle * scThetaToVFactor);
	double angle = acos(c); //using stdlib appears to be faster than cosinator :(
	double cc = cos(angle * scThetaToVFactor);

	return (0.5 + 0.5*cc);
}

#ifdef grid
double Simulation::calcInter(Spherocyl& a, Spherocyl& b) const
{
	// first, make sure they are actually close enough to interact
	const double3 p1 = a.getP();
	const double3 p2 = b.getP();

	if (space.dist2PBC(p1, p2) > scMaxInteractingCenterDistSquared) {return 0.0;}
	
	const double3 na = a.getN();
	const double3 nb = b.getN();

	if (space.spherocylDist2PBC(p1, na, p2, nb) > scMaxInteractingEndDistSquared) {return 0.0;}

	// if they are close enough to interact, check if the directions check out
	double toReturn = 0.0;

	const double3& a1 = space.getPBC(a.getS1());
	const double3& a2 = space.getPBC(a.getS2());
	const double3& b1 = space.getPBC(b.getS1());
	const double3& b2 = space.getPBC(b.getS2());

	if (scMaxTwoInters)
	{
		// only counts correctly if there can be at most one bond between two spherocyl and one sphere
		if (space.dist2PBC(a1, b1) < scMaxInteractingEndDistSquared)
		{
			double3 a1b1 = normalize(space.dispPBC(a1,b1));
			double c1 = na*a1b1;
			double c2 = -nb*a1b1;

			if (!scAngleLimited || (c1 > scCosTheta && c2 > scCosTheta)) {toReturn -= calcInterTerm(c1) * calcInterTerm(c2);}

			if (space.dist2PBC(a2, b2) < scMaxInteractingEndDistSquared)
			{
				double3 a2b2 = normalize(space.dispPBC(a2,b2));
				double d1 = -na*a2b2;
				double d2 = nb*a2b2;
				if (!scAngleLimited || (d1 > scCosTheta && d2 > scCosTheta)) {toReturn -= calcInterTerm(d1) * calcInterTerm(d2);}
			}
		}
		else if (space.dist2PBC(a1, b2) < scMaxInteractingEndDistSquared)
		{
			double3 a1b2 = normalize(space.dispPBC(a1,b2));
			double c1 = na*a1b2;
			double c2 = nb*a1b2;
			if (!scAngleLimited || (c1 > scCosTheta && c2 > scCosTheta)) {toReturn -= calcInterTerm(c1) * calcInterTerm(c2);}

			if (space.dist2PBC(a2, b1) < scMaxInteractingEndDistSquared)
			{
				double3 a2b1 = normalize(space.dispPBC(a2,b1));
				double d1 = -na*a2b1;
				double d2 = -nb*a2b1;
				if (!scAngleLimited || (d1 > scCosTheta && d2 > scCosTheta)) {toReturn -= calcInterTerm(d1)*calcInterTerm(d2);}
			}
		}
		else if (space.dist2PBC(a2, b1) < scMaxInteractingEndDistSquared)
		{
			double3 a2b1 = normalize(space.dispPBC(a2,b1));
			double c1 = -na*a2b1;
			double c2 = -nb*a2b1;
			if (!scAngleLimited || (c1 > scCosTheta && c2 > scCosTheta)) {toReturn -= calcInterTerm(c1)*calcInterTerm(c2);}
		}
		else if (space.dist2PBC(a2, b2) < scMaxInteractingEndDistSquared)
		{
			double3 a2b2 = normalize(space.dispPBC(a2,b2));
			double c1 = -na*a2b2;
			double c2 = nb*a2b2;
			if (!scAngleLimited || (c1 > scCosTheta && c2 > scCosTheta)) {toReturn -= calcInterTerm(c1)*calcInterTerm(c2);}
		}
	}
	else
	{
		if (space.dist2PBC(a1, b1) < scMaxInteractingEndDistSquared)
		{
			double3 a1b1 = normalize(space.dispPBC(a1,b1));
			double c1 = na*a1b1;
			double c2 = -nb*a1b1;

			if (!scAngleLimited || (c1 > scCosTheta && c2 > scCosTheta)) {toReturn -= calcInterTerm(c1) * calcInterTerm(c2);}
		}
		if (space.dist2PBC(a1, b2) < scMaxInteractingEndDistSquared)
		{
			double3 a1b2 = normalize(space.dispPBC(a1,b2));
			double c1 = na*a1b2;
			double c2 = nb*a1b2;
			if (!scAngleLimited || (c1 > scCosTheta && c2 > scCosTheta)) {toReturn -= calcInterTerm(c1) * calcInterTerm(c2);}
		}
		if (space.dist2PBC(a2, b1) < scMaxInteractingEndDistSquared)
		{
			double3 a2b1 = normalize(space.dispPBC(a2,b1));
			double c1 = -na*a2b1;
			double c2 = -nb*a2b1;
			if (!scAngleLimited || (c1 > scCosTheta && c2 > scCosTheta)) {toReturn -= calcInterTerm(c1)*calcInterTerm(c2);}
		}
		if (space.dist2PBC(a2, b2) < scMaxInteractingEndDistSquared)
		{
			double3 a2b2 = normalize(space.dispPBC(a2,b2));
			double c1 = -na*a2b2;
			double c2 = nb*a2b2;
			if (!scAngleLimited || (c1 > scCosTheta && c2 > scCosTheta)) {toReturn -= calcInterTerm(c1)*calcInterTerm(c2);}
		}
	}
	
	return toReturn;
	
}
#else
double Simulation::calcInter(int id1, int id2)
{
	Spherocyl& a = spherocyl[id1/2];
	Spherocyl& b = spherocyl[id2/2];

	const double3 na = a.getN();
	const double3 nb = b.getN();

	const double3& sa = space.getPBC( id1%2 == 0 ? a.getS1() : a.getS2());
	const double3& sb = space.getPBC( id2%2 == 0 ? b.getS1() : b.getS2());

	if (space.dist2PBC(sa, sb) < scMaxInteractingEndDistSquared)
	{
		double3 sasb = normalize(space.dispPBC(sa, sb));
		double c1 = (id1%2==0 ? 1.0:-1.0) * na*sasb;
		double c2 = (id2%2==0 ? -1.0:1.0) * nb*sasb;

		if (!scAngleLimited || (c1 > scCosTheta && c2 > scCosTheta))
		{
			return  -calcInterTerm(c1) * calcInterTerm(c2);
		}
	}

	return 0.0;

}
#endif

bool Simulation::interact(Spherocyl& a, Spherocyl& b)
{
	// first, make sure they are actually close enough to interact
	const double3 p1 = a.getP();
	const double3 p2 = b.getP();

	if (space.dist2PBC(p1, p2) > scMaxInteractingCenterDistSquared) {return false;}

	const double3 na = a.getN();
	const double3 nb = b.getN();

	if (space.spherocylDist2PBC(p1, na, p2, nb) > scMaxInteractingEndDistSquared) {return false;}

	// if they are close enough to interact, check if the directions check out

	const double3& a1 = space.getPBC(a.getS1());
	const double3& a2 = space.getPBC(a.getS2());
	const double3& b1 = space.getPBC(b.getS1());
	const double3& b2 = space.getPBC(b.getS2());

	if (scMaxTwoInters)
	{
		// only counts correctly if there can be at most one bond between two spherocyl and one sphere
		if (space.dist2PBC(a1, b1) < scMaxInteractingEndDistSquared)
		{
			double3 a1b1 = normalize(space.dispPBC(a1,b1));
			double c1 = na*a1b1;
			double c2 = -nb*a1b1;

			if (!scAngleLimited || (c1 > scCosTheta && c2 > scCosTheta)) {return true;}

			if (space.dist2PBC(a2, b2) < scMaxInteractingEndDistSquared)
			{
				double3 a2b2 = normalize(space.dispPBC(a2,b2));
				double d1 = -na*a2b2;
				double d2 = nb*a2b2;
				if (!scAngleLimited || (d1 > scCosTheta && d2 > scCosTheta)) {return true;}
			}
		}
		else if (space.dist2PBC(a1, b2) < scMaxInteractingEndDistSquared)
		{
			double3 a1b2 = normalize(space.dispPBC(a1,b2));
			double c1 = na*a1b2;
			double c2 = nb*a1b2;
			if (!scAngleLimited || (c1 > scCosTheta && c2 > scCosTheta)) {return true;}

			if (space.dist2PBC(a2, b1) < scMaxInteractingEndDistSquared)
			{
				double3 a2b1 = normalize(space.dispPBC(a2,b1));
				double d1 = -na*a2b1;
				double d2 = -nb*a2b1;
				if (!scAngleLimited || (d1 > scCosTheta && d2 > scCosTheta)) {return true;}
			}
		}
		else if (space.dist2PBC(a2, b1) < scMaxInteractingEndDistSquared)
		{
			double3 a2b1 = normalize(space.dispPBC(a2,b1));
			double c1 = -na*a2b1;
			double c2 = -nb*a2b1;
			if (!scAngleLimited || (c1 > scCosTheta && c2 > scCosTheta)) {return true;}
		}
		else if (space.dist2PBC(a2, b2) < scMaxInteractingEndDistSquared)
		{
			double3 a2b2 = normalize(space.dispPBC(a2,b2));
			double c1 = -na*a2b2;
			double c2 = nb*a2b2;
			if (!scAngleLimited || (c1 > scCosTheta && c2 > scCosTheta)) {return true;}
		}
	}
	else
	{
		if (space.dist2PBC(a1, b1) < scMaxInteractingEndDistSquared)
		{
			double3 a1b1 = normalize(space.dispPBC(a1,b1));
			double c1 = na*a1b1;
			double c2 = -nb*a1b1;

			if (!scAngleLimited || (c1 > scCosTheta && c2 > scCosTheta)) {return true;}
		}
		if (space.dist2PBC(a1, b2) < scMaxInteractingEndDistSquared)
		{
			double3 a1b2 = normalize(space.dispPBC(a1,b2));
			double c1 = na*a1b2;
			double c2 = nb*a1b2;
			if (!scAngleLimited || (c1 > scCosTheta && c2 > scCosTheta)) {return true;}
		}
		if (space.dist2PBC(a2, b1) < scMaxInteractingEndDistSquared)
		{
			double3 a2b1 = normalize(space.dispPBC(a2,b1));
			double c1 = -na*a2b1;
			double c2 = -nb*a2b1;
			if (!scAngleLimited || (c1 > scCosTheta && c2 > scCosTheta)) {return true;}
		}
		if (space.dist2PBC(a2, b2) < scMaxInteractingEndDistSquared)
		{
			double3 a2b2 = normalize(space.dispPBC(a2,b2));
			double c1 = -na*a2b2;
			double c2 = nb*a2b2;
			if (!scAngleLimited || (c1 > scCosTheta && c2 > scCosTheta)) {return true;}
		}
	}

	return false;

}

StaticCluster Simulation::getCluster(Spherocyl& s)
{
	StaticCluster toReturn;
	toReturn.push_back(s.getID());

	vector<bool> isAccountedFor(scN, false);
	isAccountedFor[s.getID()] = true;

	vector<int> toCheckThisRound;
	vector<int> toCheckNextRound;

	toCheckThisRound.push_back(s.getID());

	while (toCheckThisRound.size() != 0)
	{
		for (const auto& startID : toCheckThisRound)
		{
#ifdef grid
			space.getIDsNear(spherocyl[startID].getP());

			while (space.gotMoreIDs())
			{
				int possiblePartnerID = space.getNextID();
				if (!isAccountedFor[possiblePartnerID] && interact(spherocyl[startID], spherocyl[possiblePartnerID]))
				{
					toCheckNextRound.push_back(possiblePartnerID);
					isAccountedFor[possiblePartnerID] = true;
					toReturn.push_back(possiblePartnerID);
				}
			}
#else
			// Kd-tree
			for (int side = 0; side < 2; side++)
			{

				int tipID = 2*startID + side;
				auto box = getInteractionBox(tipID);
				auto ids = space.getInBoxTip(box.first, box.second);

				for (const auto& possiblePartnerTipID : ids)
				{
					int possiblePartnerID = possiblePartnerTipID/2;
					if (!isAccountedFor[possiblePartnerID] && calcInter(tipID, possiblePartnerTipID) < 0.0)
					{
						toCheckNextRound.push_back(possiblePartnerID);
						isAccountedFor[possiblePartnerID] = true;
						toReturn.push_back(possiblePartnerID);
					}
				}
			}
#endif
		}
		toCheckThisRound = toCheckNextRound;
		toCheckNextRound.clear();
	}

	return toReturn;
}


void Simulation::calcE()
{
	// reset E
	E = 0.0;

	// calculate all interaction terms
#ifdef grid
	for (int id1 = 0; id1 < scN; id1++)
	{
		auto& s = spherocyl[id1];

		space.getIDsNear(s.getP());
		while (space.gotMoreIDs())
		{
			int id2 = space.getNextID();

			if (id1 < id2) // avoid double-counting
			{
				E += calcInter(s, spherocyl[id2]);
			}
		}
	}
#else
	// Kd-tree
	for (int idS = 0; idS < 2*scN; idS++)
	{
		auto box = getInteractionBox(idS);
		auto ids = space.getInBoxTip(box.first, box.second);

		for (const auto& id2 : ids)
		{
			if (idS < id2) // avoid double-counting
			{
				E += calcInter(idS, id2);
			}
		}
	}
#endif
}

double Simulation::calcEnergyContribution(Spherocyl& s)
{
	double toReturn = 0.0;

	int id1 = s.getID();
#ifdef grid
	space.getIDsNear(s.getP());

	while (space.gotMoreIDs())
	{
		int id2 = space.getNextID();

		if (id1 != id2)
		{
			toReturn += calcInter(s, spherocyl[id2]);
		}
	}
#else
	// Kd-tree: calculate energy at both end points
	for (int i = 0; i < 2; i++)
	{
		int idS = 2 * s.getID() + i;
		auto box = getInteractionBox(idS);
		auto ids = space.getInBoxTip(box.first, box.second);
		for (const auto& id2 : ids)
		{
			if (idS != id2)
			{
				toReturn += calcInter(idS, id2);
			}
		}
	}
#endif

	return toReturn;
}
bool Simulation::interactsWithSomeone(Spherocyl& s)
{
	int id1 = s.getID();
#ifdef grid
	space.getIDsNear(s.getP());

	while (space.gotMoreIDs())
	{
		int id2 = space.getNextID();

		if (id1 != id2)
		{
			if (interact(s, spherocyl[id2]))
			{
				return true;
			}
		}
	}
	return false;
#else
	auto box = getIntersectionBox(s);
	return !space.traverse(&Simulation::interactsWithSomeoneTraverse, box.first, box.second);
#endif
}

void Simulation::setState(const SystemState<GridSize, Spherocyl, SpherocylParameter>& state)
{
	// load parameters
	calcSpherocylParameters(state.extra);

	// load thermodynamical properties
	// E = state.E; // not required since we recalculate E
	beta = state.beta;
#ifdef NPT
	P = state.P;
#endif

	// set up grid
#ifdef grid
	space.initialize(state.size.gridSize, scMinGridSize, state.size.spaceSize);
#else
	space.initialize(state.size, this);
#endif
	// copy spherocylinders incl positions
	spherocyl = state.particle;
	scN = spherocyl.size();

	// make sure everything is within PBC and normalized
	for (int i = 0; i < scN; i++)
	{
		Spherocyl& s = spherocyl[i];
		s.setP(space.getPBC(s.getP()));
		s.fixRoundingErrors();
	}

	// register everyting on grid
	for (int i = 0; i < scN; i++)
	{
		space.insert(spherocyl[i]);
	}

	// set system energy
	calcE();
}

SystemState<GridSize, Spherocyl, SpherocylParameter> Simulation::getState()
{
	SystemState<GridSize, Spherocyl, SpherocylParameter> state;

	state.particle = spherocyl;
#ifdef grid
	state.size.gridSize = space.getMaxCells();
	state.size.spaceSize = space.getSpaceSize();
#else
	state.size = space.getSpaceSize();
#endif

	state.E = E;
	state.beta = beta;
#ifdef NPT
	state.P = P;
#endif

	state.extra = spherocylParameters;

	return state;
}

double Simulation::getE() const
{
	return E;
}

double Simulation::getV() const
{
	double3 s = space.getSpaceSize();
	return s.x*s.y*s.z;
}

void Simulation::setBeta(const double& b)
{
	beta = b;
}

#ifdef NPT
void Simulation::setPressure(const double& p)
{
	P = p;
}
#endif

bool Simulation::decide(const double& dE)
{
	if (dE <= 0.0 || beta == 0.0) {return true;}

	return (random->random() < exp(-beta*dE));
}

bool Simulation::MCSingleGlobalTranslationRotation(const void* parameter)
{
	// pick a spherocylinder
	auto& s = spherocyl[random->random(spherocyl.size())];

	// propose transition
	double3 size = space.getSpaceSize();
	double3 newPos = {random->random()*size.x, random->random()*size.y, random->random()*size.z};
	double3 newN = random->randomSpherical();
	double3 oldPos = s.getP();
	double3 oldN = s.getN();

	return MCSinglePN(s, newPos, newN);
}
bool Simulation::MCSingleGlobalTranslation(const void* parameter)
{
	// pick a spherocylinder
	auto& s = spherocyl[random->random(spherocyl.size())];

	// propose transition
	double3 size = space.getSpaceSize();
	double3 newPos = {random->random()*size.x, random->random()*size.y, random->random()*size.z};
	double3 oldPos = s.getP();


	return MCSinglePN(s, newPos, s.getN());
}

bool Simulation::MCSingleLocalTranslation(const void* parameter)
{
	// pick a spherocylinder
	auto& s = spherocyl[random->random(spherocyl.size())];

	// propose transition
	double3 oldPos = s.getP();
	double3 newPos = space.getPBC( oldPos + 0.1*random->randomSpherical() );

	return MCSinglePN(s, newPos, s.getN());
}

bool Simulation::MCSingleRotation(const void* parameter)
{
	// pick a spherocylinder
	auto& s = spherocyl[random->random(spherocyl.size())];

	// propose transition
	double3 newN = random->randomSpherical();
	double3 oldPos = s.getP();
	double3 oldN = s.getN();

	return MCSinglePN(s, oldPos, newN);
}

bool Simulation::MCSingleLocalRotation(const void* parameter)
{
	// pick a spherocylinder
	auto& s = spherocyl[random->random(spherocyl.size())];

	// propose transition
	double3 oldPos = s.getP();
	double3 oldN = s.getN();

	bool which = (random->random() < 0.5);
	double angle = (random->random()*2 -1.0)*30;

	double3 origin = (which ? s.getS1() : s.getS2());
	double3 axis = random -> randomSpherical();
	double3 newN = fp::rotate(oldN, axis, angle);
	double3 s1disp = newN - oldN;
	double3 newPos = (which ? space.getPBC(oldPos - s1disp) : space.getPBC(oldPos + s1disp) );

	return MCSinglePN(s, newPos, newN);
}

bool Simulation::MCSingleNormalTranslation(const void* parameter)
{
	// pick a spherocylinder
	auto& s = spherocyl[random->random(spherocyl.size())];

	// propose transition
	double3 oldPos = s.getP();
	double3 newPos = space.getPBC( oldPos + 2.0*(random->random()-0.5)*spherocylParameters.d*s.getN() );

	return MCSinglePN(s, newPos, s.getN());
}

bool Simulation::MCSinglePN(Spherocyl& s, const double3& newPos, const double3& newN)
{
	double3 oldPos = s.getP();
	double3 oldN = s.getN();

	// check if place available
	space.remove(s);
	s.setPN(newPos, newN);

	if (intersects(s))
	{
		s.setPN(oldPos, oldN);
		space.insert(s);
		return false;
	}

	// calculate energy difference
	double Eafter = calcEnergyContribution(s);
	s.setPN(oldPos, oldN);
	double Ebefore = calcEnergyContribution(s);

	if (decide(Eafter - Ebefore))
	{
		// accept
		s.setPN(newPos, newN);
		space.insert(s);
		E+=Eafter-Ebefore;
		return true;
	}
	else
	{
		// reject: nothing to do
		space.insert(s);
		return false;
	}

}

#ifdef NPT
bool Simulation::MCVolumeIso(const void* parameter)
{
	double3 oldSpaceSize = space.getSpaceSize();

	double V0 = oldSpaceSize.x*oldSpaceSize.y*oldSpaceSize.z;
	double volumeScaleFactor = 0.995 + 0.01*random->random();
	double V1 = V0 * volumeScaleFactor;
	double linearScaleFactor = std::pow(volumeScaleFactor, (1.0/3.0) );
	double3 newSpaceSize = oldSpaceSize*linearScaleFactor;
	if (newSpaceSize.x <= scMinSpaceDimension || newSpaceSize.y <= scMinSpaceDimension || newSpaceSize.z <= scMinSpaceDimension)
	{
		return false;
	}

	// calculate energy difference
	double oldE = E;

	// save old positions
	vector<double3> oldP(scN);
	vector<double3> newP(scN);
	for (int i = 0; i < scN; i++)
	{
		oldP[i] = spherocyl[i].getP();
	}


	// lift everything from grid
#ifdef grid
	for (int i = 0; i < scN; i++)
	{
		space.remove(spherocyl[i]);
	}
#else
	space.clear();
#endif

	// resize grid
	space.changeVolume(newSpaceSize);


	// place on new positions
	bool worksSoFar = true;
	for (int i = 0; i < scN; i++)
	{
		newP[i] = oldP[i] * linearScaleFactor;
		spherocyl[i].setP( newP[i] );
		if (intersects(spherocyl[i]))
		{
			// already occupied. abort
			worksSoFar = false;

			for (int j = 0; j < i; j++)
			{
				space.remove(spherocyl[j]);
			}

			break;
		}
		else
		{
			space.insert(spherocyl[i]);
		}
	}

	if (!worksSoFar) // space too occupied
	{
		space.changeVolume(oldSpaceSize);

		for (int i = 0; i < scN; i++)
		{
			spherocyl[i].setP(oldP[i]);
			space.insert(spherocyl[i]);
		}

		return false;
	}

	// space not too occupied. all particles are placed in new pos. calc E difference
	calcE();
	double newE = E;
	double dE = newE-oldE;
	double dV = V1-V0;

	double accP = exp( fmin( 0, (scN)*log(V1/V0)-beta*dE-beta*P*dV ) );
	if (accP >= 1.0 || random->random() < accP)
	{
		// accept
		return true;
	}
	else
	{
		// reject
#ifdef grid
		for (int i = 0; i < scN; i++)
		{
			space.remove(spherocyl[i]);
		}
#else
		space.clear();
#endif
		space.changeVolume(oldSpaceSize);
		for (int i = 0; i < scN; i++)
		{
			spherocyl[i].setP(oldP[i]);
			space.insert(spherocyl[i]);
		}
		E = oldE;
		return false;
	}
}
#endif

#ifdef NPT
bool Simulation::MCVolumeAniso(const void* parameter)
{
	double3 oldSpaceSize = space.getSpaceSize();

	double V0 = oldSpaceSize.x*oldSpaceSize.y*oldSpaceSize.z;
	double linearScaleFactor = 0.995 + 0.01*random->random(); //TOD: use parameter
	double V1 = V0 * linearScaleFactor;
	double3 newSpaceSize = oldSpaceSize;

	int dir = random->random(3);
	switch (dir)
	{
		case 0:
			newSpaceSize.x *= linearScaleFactor;
			break;
		case 1:
			newSpaceSize.y *= linearScaleFactor;
			break;
		case 2:
			newSpaceSize.z *= linearScaleFactor;
			break;
	}
	if (newSpaceSize.x <= scMinSpaceDimension || newSpaceSize.y <= scMinSpaceDimension || newSpaceSize.z <= scMinSpaceDimension)
	{
		return false;
	}

	// check aspect ratio
	const double maxAR = 2.0;
	if (newSpaceSize.x / newSpaceSize.y > maxAR || newSpaceSize.x / newSpaceSize.y < 1.0/maxAR) {return false;}
	if (newSpaceSize.y / newSpaceSize.z > maxAR || newSpaceSize.y / newSpaceSize.z < 1.0/maxAR) {return false;}
	if (newSpaceSize.z / newSpaceSize.x > maxAR || newSpaceSize.z / newSpaceSize.x < 1.0/maxAR) {return false;}

	// calculate energy difference
	double oldE = E;

	// save old positions
	vector<double3> oldP(scN);
	vector<double3> newP(scN);
	for (int i = 0; i < scN; i++)
	{
		oldP[i] = spherocyl[i].getP();
	}


	// lift everything from grid
#ifdef grid
	for (int i = 0; i < scN; i++)
	{
		space.remove(spherocyl[i]);
	}
#else
	space.clear();
#endif

	// resize grid
	space.changeVolume(newSpaceSize);


	// place on new positions
	bool worksSoFar = true;
	for (int i = 0; i < scN; i++)
	{
		newP[i] = oldP[i];

		switch (dir)
		{
			case 0:
				newP[i].x *= linearScaleFactor;
				break;
			case 1:
				newP[i].y *= linearScaleFactor;
				break;
			case 2:
				newP[i].z *= linearScaleFactor;
				break;
		}


		spherocyl[i].setP( newP[i] );
		if (intersects(spherocyl[i]))
		{
			// already occupied. abort
			worksSoFar = false;

			for (int j = 0; j < i; j++)
			{
				space.remove(spherocyl[j]);
			}

			break;
		}
		else
		{
			space.insert(spherocyl[i]);
		}
	}

	if (!worksSoFar) // space too occupied
	{
		space.changeVolume(oldSpaceSize);

		for (int i = 0; i < scN; i++)
		{
			spherocyl[i].setP(oldP[i]);
			space.insert(spherocyl[i]);
		}

		return false;
	}

	// space not too occupied. all particles are placed in new pos. calc E difference
	calcE();
	double newE = E;
	double dE = newE-oldE;
	double dV = V1-V0;

	double accP = exp( fmin( 0, (scN)*log(V1/V0)-beta*dE-beta*P*dV ) );
	if (accP >= 1.0 || random->random() < accP)
	{
		// accept
		return true;
	}
	else
	{
		// reject
#ifdef grid
		for (int i = 0; i < scN; i++)
		{
			space.remove(spherocyl[i]);
		}
#else
		space.clear();
#endif
		space.changeVolume(oldSpaceSize);
		for (int i = 0; i < scN; i++)
		{
			spherocyl[i].setP(oldP[i]);
			space.insert(spherocyl[i]);
		}
		E = oldE;
		return false;
	}
}
#endif

bool Simulation::MCClusterMove(const void*)
{
	// pick a spherocylinder
	auto& s = spherocyl[random->random(spherocyl.size())];

	auto clusterBefore = getCluster(s);
	if (clusterBefore.size() == scN) {return false;} // shifting all spherocylinders makes no sense


	double3 shift = 0.1*random->randomSpherical();

	// save old stuff
	vector<double3> oldPos(clusterBefore.size());
	vector<double3> newPos(clusterBefore.size());
	for (int i = 0; i < clusterBefore.size(); i++)
	{
		oldPos[i] = spherocyl[clusterBefore[i]].getP();
	}

	// check if there is room
	for (int i = 0; i < clusterBefore.size(); i++)
	{
		int id = clusterBefore[i];
		space.remove(spherocyl[id]);

	}

	for (int i = 0; i < clusterBefore.size(); i++)
	{
		int id = clusterBefore[i];
		newPos[i] = space.getPBC(oldPos[i] + shift);
		spherocyl[id].setP(newPos[i]);

		if (intersects(spherocyl[id]))
		{
			// restore all the things
			for (int j = 0; j <= i; j++)
			{
				spherocyl[clusterBefore[j]].setP(oldPos[j]);
			}
			for (int j = 0; j < clusterBefore.size(); j++)
			{
				space.insert(spherocyl[clusterBefore[j]]);
			}

			return false;
		}
	}

	// looks like there is room. see if it would be interacting (this changes cluster size, so abort!)
	for (const auto& id : clusterBefore)
	{
		if (calcEnergyContribution(spherocyl[id]) != 0.0)
		{
			// reject!
			for (int i = 0; i < clusterBefore.size(); i++)
			{
				spherocyl[clusterBefore[i]].setP(oldPos[i]);
				space.insert(spherocyl[clusterBefore[i]]);
			}
			return false;
		}
	}

	// accept
	for (int i = 0; i < clusterBefore.size(); i++)
	{
		space.insert(spherocyl[clusterBefore[i]]);
	}

	return true;
}



void Simulation::checkSpace()
{
#ifdef grid
	// build second grid
	Grid g2;
	double r,d,t;
	r = spherocylParameters.r;
	d = spherocylParameters.d;
	t = spherocylParameters.t;

	g2.initialize(space.getMaxCells(), 2.0+2*r+d, space.getSpaceSize());

	for (int i = 0; i < scN; i++)
	{
		g2.insert(spherocyl[i]);
	}

	// compare to given grid
	for (int i = 0; i < scN; i++)
	{
		space.getIDsNear(spherocyl[i].getP());
		g2.getIDsNear(spherocyl[i].getP());

		while (space.gotMoreIDs())
		{
			if (!g2.gotMoreIDs())
			{
				cerr << "diff length!" << endl;
			}
			int ida = space.getNextID();
			int idb = g2.getNextID();
			if (ida != idb)
			{
				cerr << "not same id near " << spherocyl[i].getP() << ": " << ida << " vs " << idb	<< endl;
			}
		}
	}

	double fakeE = E;
	calcE();
	double realE = E;

	if (realE != fakeE)
	{
		cerr << realE << " != " << fakeE << endl;
	}
#else
	KdTreeSpherocyl g2;
	double r,d,t;
	r = spherocylParameters.r;
	d = spherocylParameters.d;
	t = spherocylParameters.t;

	g2.initialize(space.getSpaceSize(), this);

	for (int i = 0; i < scN; i++)
	{
		g2.insert(spherocyl[i]);
	}

	// compare to given grid
	for (int i = 0; i < scN; i++)
	{
		auto box = getIntersectionBox(spherocyl[i]);
		/*cout << spherocyl[i].getP() << endl;
		  cout << spherocyl[i].getN() << endl;
		  cout << box.first << endl;
		  cout << box.second << endl;*/


		auto ids1 = space.getInBoxCOM(box.first, box.second);
		auto ids2 = g2.getInBoxCOM(box.first, box.second);

		if (ids2.size() == 0)
		{
			cout << spherocyl[i].getP() << endl << box.first << endl << box.second << endl;
		}
		if (ids1.size() != ids2.size())
		{
			cerr << "there should be " << ids2.size() << ", but there are " << ids1.size() << endl;
		}
		auto ids3 = space.getInBoxTip(box.first, box.second);
		auto ids4 = g2.getInBoxTip(box.first, box.second);
		if (ids4.size() == 0)
		{
			cout << spherocyl[i].getP() << endl << box.first << endl << box.second << endl;
		}
		if (ids3.size() != ids4.size())
		{
			cerr << "there should be " << ids4.size() << ", but there are " << ids3.size() << endl;
		}
	}

	double fakeE = E;
	calcE();
	double realE = E;

	if (realE != fakeE)
	{
		cerr << realE << " != " << fakeE << endl;
	}

#endif

}

#ifndef grid
bool Simulation::intersectsTraverse(int id2)
{
	return !(globalId1 != id2 && intersect(spherocyl[globalId1], spherocyl[id2]));
}
bool Simulation::interactsWithSomeoneTraverse(int id2)
{
	return !(globalId1 != id2 && interact(spherocyl[globalId1], spherocyl[id2]));
}

pair<double3, double3> Simulation::getIntersectionBox(Spherocyl& s)
{
	double3 l, r, t1, t2;

	t1 = s.getP() + s.getN();
	t2 = s.getP() - s.getN();

	l.x = fmin(t1.x, t2.x) - intersectionBoxBorderSize;
	l.y = fmin(t1.y, t2.y) - intersectionBoxBorderSize;
	l.z = fmin(t1.z, t2.z) - intersectionBoxBorderSize;

	r.x = fmax(t1.x, t2.x) + intersectionBoxBorderSize;
	r.y = fmax(t1.y, t2.y) + intersectionBoxBorderSize;
	r.z = fmax(t1.z, t2.z) + intersectionBoxBorderSize;

	pair<double3, double3> toReturn;
	toReturn.first = space.getPBC(l);
	toReturn.second = space.getPBC(r);

	return toReturn;
}

pair<double3, double3> Simulation::getInteractionBox(int idS)
{
	double3 l, r, s;

	s = ( idS%2==0 ? spherocyl[idS/2].getS1() : spherocyl[idS/2].getS2() );


	l.x = s.x - interactionBoxBorderSize;
	l.y = s.y - interactionBoxBorderSize;
	l.z = s.z - interactionBoxBorderSize;

	r.x = s.x + interactionBoxBorderSize;
	r.y = s.y + interactionBoxBorderSize;
	r.z = s.z + interactionBoxBorderSize;

	pair<double3, double3> toReturn;
	toReturn.first = space.getPBC(l);
	toReturn.second = space.getPBC(r);

	return toReturn;
}
#endif
