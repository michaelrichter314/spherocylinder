#include "icmaker.h"

ICmaker::ICmaker(Randomizer* R) : r(R)
{
	type = ICrandom;
}

ICmaker::ICmaker()
{
	type = ICsquarelattice;
}

ICmaker::~ICmaker()
{
	if (itsMyRandomizer)
	{
		delete r;
	}
}

void ICmaker::setGridSize(const int3& gs)
{
	gridSize = gs;
}

void ICmaker::setSpaceSize(const double3& ss)
{
	spaceSize = ss;
}

void ICmaker::setSpherocylParameters(const SpherocylParameter& SP)
{
	sp = SP;
}

void ICmaker::setN(unsigned int n)
{
	N = n;
}

void ICmaker::setBeta(const double& b)
{
	beta = b;
}

#ifdef NPT
void ICmaker::setP(const double& p)
{
	P = p;
}
#endif

void ICmaker::setType(const ICtype& t)
{
	type = t;
}

SystemState<GridSize, Spherocyl, SpherocylParameter> ICmaker::getState()
{
	make();
	return state;
}

void ICmaker::make()
{
	if (hasBeenMade) {return;}

	// initialize grid
#ifdef grid
	space.initialize(gridSize, 2.0+2*sp.r + sp.d, spaceSize);
#else
	space.initialize(spaceSize, NULL);
#endif
	
	// make sure we have a randomizer
	if (type == ICrandom && r == NULL)
	{
		itsMyRandomizer = true;
		r = new Randomizer();

		cerr << "warning: ICmaker has not been given a randomizer\n";
	}

	// make room for spherocylinders in memory
	spherocyl.resize(N);

	// set pherocyl ID's
	for (int i = 0; i < N; i++)
	{
		spherocyl[i].setID(i);
	}

	// place on grid
	if (type == ICrandom)
	{
		unsigned long trials = 0;
		unsigned long maxTrials = N * 100'000; //TODO: const number somewhere
		
		do
		{
#ifdef grid
			space.initialize(gridSize, 2.0+2*sp.r + sp.d, spaceSize);
#else
			space.initialize(spaceSize, NULL);
#endif
			for (int i = 0; i < N; i++)
			{
				cerr << "placing " << i << endl;
				double3 newP;
				double3 newN;

				do
				{
					newP = { spaceSize.x*r->random(), spaceSize.y*r->random(), spaceSize.z*r->random() };
					newN = r -> randomSpherical();

					spherocyl[i].setPN(newP, newN);

					trials++;
					if (trials >= maxTrials) {break;}
				} while (intersectsSomeone(spherocyl[i]));

				if (trials >= maxTrials)
				{
					spherocyl.clear();
					break;
				}

				space.insert(spherocyl[i]);
			}
		} while (spherocyl.size() != N);

	}
	else
	{
		//TODO
		cerr << "error: non-radom inicital condition not supported\n";
	}

	// make state
	state.particle = spherocyl;
#ifdef grid
	state.size = {gridSize, spaceSize};
#else
	state.size = spaceSize;
#endif
	state.E = 0.0;
	state.beta = beta;
#ifdef NPT
	state.P = P;
#endif
	state.extra = sp;
}

bool ICmaker::intersectsSomeone(Spherocyl& s)
{
	int id1 = s.getID();
#ifdef grid
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
	auto box = getIntersectionBox(s);
	auto ids = space.getInBoxCOM(box.first, box.second);

	for (const auto& id2 : ids)
	{
		if (id1 != id2 && intersect(s, spherocyl[id2]))
		{
			return true;
		}

	}

	return false;
#endif
	
}


bool ICmaker::intersect(Spherocyl& a, Spherocyl& b) const
{
	double3 p1 = a.getP();
	double3 p2 = b.getP();

	// if they are close, calculate actual distance
	if (space.dist2PBC(p1, p2) <= (2*sp.r+2)*(2*sp.r+2))
	{
		// check actual distance
		if (space.spherocylDist2PBC(p1, a.getN(), p2, b.getN()) <= 4*sp.r*sp.r)
		{
			return true;
		}
	}
	return false;
}

#ifndef grid
pair<double3, double3> ICmaker::getIntersectionBox(Spherocyl& s)
{
	const double intersectionBoxBorderSize = 1.0 + sp.r + sp.d;
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
#endif
	

