#include "analyzer.h"

Analyzer::Analyzer(string filename)
{
	loader = make_unique<VTFloader>(filename);
	timesteps = loader -> getHandler() -> getTimesteps();
	timestep = 0;
	setState();
}

void Analyzer::setTheta(double t)
{
	theta = t;
}

bool Analyzer::goToNextTimestep()
{
	timestep++;
	setState();
	return (timestep <= timesteps);
}

double Analyzer::getOrderParameter()
{
	double sum = 0;
	int n = 0;

	for (int i = 0; i < state.particle.size(); i++)
	{
		for (int j = i+1; j < state.particle.size(); j++)
		{
			double c = state.particle[i].getN() * state.particle[j].getN();
			sum += (3*c*c-1)/2.0;
			n++;
		}
	}

	return sum/n;
}

double Analyzer::getSmecticParameter()
{
	vector<double3> trirector;

	for (int i = 0; i < state.particle.size(); i++)
	{
		trirector.push_back(getTrirector(state.particle[i]));
		//std::cout << getTrirector(state.particle[i]) << std::endl;
	}

	double sum = 0;
	int n = 0;

	vector<double> cosss;
	for (int i = 0; i < trirector.size(); i++)
	{
		for (int j = i+1; j < trirector.size(); j++)
		{
			double c = trirector[i]*trirector[j];
			sum += (3*c*c-1)/2.0;
			cosss.push_back(c*c);
			//std::cout << c << " " << trirector[i] << " " << trirector[j] << " " << sum << std::endl;
			n++;
		}
	}

	//return variance(&cosss[0], cosss.size());

	return sum/n;
}

double Analyzer::getAvNN()
{
	double NN = 0.0;
	for (int i = 0; i < state.particle.size(); i++)
	{
		//double d = getNN(state.particle[i]);
		//std::cout << i << " " << d << std::endl;
		NN += getNN(state.particle[i]);
	}

	return NN / state.particle.size();
}

double Analyzer::getSmect()
{
	int NN = 0;
	for (int i = 0; i < state.particle.size(); i++)
	{
		//double d = getNN(state.particle[i]);
		//std::cout << i << " " << d << std::endl;
		NN += countSmect(state.particle[i]);
	}

	return (double)NN / state.particle.size();
}

double Analyzer::getNN(Spherocyl& s) //not quite correct due to box stuff
{
	double dist = 0.5;
	while (true)
	{
		vector<int> candidates;
		if (2*dist >= state.size.x || 2*dist >= state.size.y || 2*dist > state.size.z)
		{
			// look in whole box
			for (int i = 0; i < state.particle.size(); i++)
			{
				candidates.push_back(i);
			}
		}
		else
		{
			// candidates in neighborhood

			double3 ll = s.getP();
			ll.x -= dist;
			ll.y -= dist;
			ll.z -= dist;
			double3 ur = s.getP();
			ur.x += dist;
			ur.y += dist;
			ur.z += dist;

			ll = space.getPBC(ll);
			ur = space.getPBC(ur);

			candidates = space.getInBoxCOM(ll, ur);
		}

		double lowest;
		bool foundOne = false;
	   if (s.getID() != 0)
	   {
		   lowest = space.dist2PBC(s.getP(), state.particle[0].getP());
	   }
	   else
	   {
		   lowest = space.dist2PBC(s.getP(), state.particle[1].getP());
	   }

	   for (int i = 0; i < candidates.size(); i++)
	   {
		   int id = candidates[i];
		   if (id != s.getID())
		   {
			   double d = space.dist2PBC(s.getP(), state.particle[id].getP());
			   if (d <= lowest) {lowest = d; foundOne = true;}
		   }
	   }

	   if (foundOne) {return sqrt(lowest);}

	   dist *= 2.0;
	}
}

int Analyzer::countSmect(Spherocyl& s)
{
	int id = s.getID();
	auto p = s.getP();
	auto N = s.getN();
	
	int count = 0;
	for (int i = 0; i < state.particle.size(); i++)
	{
		double d = space.dist2PBC(p, state.particle[i].getP());
		if (d < 1.0)
		{
			// calc pos along normal
			double dd = (state.particle[i].getP() - p) * N;

			if (fabs(dd) < 0.2) {count++;}
		}
	}

	return count;
}

double3 Analyzer::getTrirector(Spherocyl& s)
{
	double3 p0 = s.getP();
	double3 p1;
	double3 p2;

	double nearest = -1;
	double secondnearest = -1;

	for (int i = 0; i < state.particle.size(); i++)
	{
		if (i != s.getID())
		{
			auto& s1 = state.particle[i];

			double d = space.dist2PBC(p0, s1.getP());
			if ( nearest < 0 || d < nearest )
			{
				// shift 2nd
				secondnearest = nearest;
				p2 = p1;

				p1 = s1.getP();
				nearest = d;
			}
			else if (secondnearest < 0 || d < secondnearest)
			{
				p2 = s1.getP();
				secondnearest = d;
			}
		}
	}

	// now we got 3 points. generate normal
	double3 d1 = space.dispPBC(p0, p1);
	double3 d2 = space.dispPBC(p0, p2);

	return normalize(cross(d1, d2));
}



void Analyzer::setState()
{
	state = loader -> getState(timestep);
	space.clear();
	space.changeVolume(state.size);
	for (int i = 0; i < state.particle.size(); i++)
	{
		space.insert(state.particle[i]);
	}

	state.extra.t = theta;
}

pair<double3, double3> Analyzer::getInteractionBox(int idS)
{
	double3 l, r, s;
	double interactionBoxBorderSize = 2*state.extra.r+state.extra.d;
	s = ( idS%2==0 ? state.particle[idS/2].getS1() : state.particle[idS/2].getS2() );


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

bool Analyzer::interact(int id1, int id2)
{
	if (id1 == id2) {return false;}

	double r = state.extra.r;
	double d = state.extra.d;
	double t = state.extra.t;
	double scMaxInteractingEndDistSquared = (2+2*r+d)*(2+2*r+d);
	double scAngleLimited = ( t < 180.0 );
	const double PI = 3.14159265359;
	double scCosTheta = cos(t*PI/180.0);

	Spherocyl& a = state.particle[id1/2];
	Spherocyl& b = state.particle[id2/2];

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
			return true;
		}
	}

	return false;
}

double Analyzer::getAssociationSpherocyl()
{
	int assos = 0;

	for (auto&& s : state.particle)
	{
		int id = s.getID();
		int tipID = 2*id + 0;
		auto box = getInteractionBox(tipID);
		auto ids = space.getInBoxTip(box.first, box.second);
		bool foundOne = false;
		for (const auto& pTipID : ids)
		{
			if (interact(tipID, pTipID))
			{
				assos++;
				foundOne = true;
				break;
			}
		}

		if (!foundOne)
		{
			tipID = 2*id + 1;
			box = getInteractionBox(tipID);
			ids = space.getInBoxTip(box.first, box.second);
			for (const auto& pTipID : ids)
			{
				if (interact(tipID, pTipID))
				{
					assos++;
					break;
				}
			}
		}
	}

	return (double)assos / state.particle.size();
}
double Analyzer::getAssociationEnd()
{
	int assos = 0;

	for (auto&& s : state.particle)
	{
		int id = s.getID();
		int tipID = 2*id + 0;
		auto box = getInteractionBox(tipID);
		auto ids = space.getInBoxTip(box.first, box.second);
		for (const auto& pTipID : ids)
		{
			if (interact(tipID, pTipID))
			{
				assos++;
				break;
			}
		}

		tipID = 2*id + 1;
		box = getInteractionBox(tipID);
		ids = space.getInBoxTip(box.first, box.second);
		for (const auto& pTipID : ids)
		{
			if (interact(tipID, pTipID))
			{
				assos++;
				break;
			}
		}
	}

	return (double)assos / (2*state.particle.size());
}
