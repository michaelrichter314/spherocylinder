#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fp/fpdim3.h>
#include <fp/fprandom.h>
#include <fp/fpStatistics.h>

#define NOFEEDBACK
//#define SPOTS

const double PI = 3.14159265359;

using namespace std;
using namespace fp;
//using namespace lcmc;

const double MCmoveDist = 0.28;
const double MCrotAngle = 4.0;
const double MCcluMoveDist = 1.0;
const double MCcluRotAngle = 6.0;
#ifdef SPOTS
	const double MCspotDist = 0.1;
#endif

struct SimSettings
{
	unsigned int sc;
	double scR; // radius of spherocylindder
	double scD; // distance of interaction

#ifdef SPOTS
	int scSpots;
	double scK;
#endif

	double scOneOverT; //anglge theta

	int3 size; // as multiple of gridSize
};

#ifdef NOFEEDBACK
typedef void MCfeedback;
typedef void IterationFeedback;
#else
struct MCfeedback
{
	double dist = 0.0;
};
struct IterationFeedback
{
	double move = 0.0;
	double rot = 0.0;
	double cluMove = 0.0;
	double cluRot = 0.0;
#ifdef SPOTS
	double spotMove = 0.0;
#endif
};
#endif



class Simulation
{
	protected:
		struct Spherocyl
		{
			Spherocyl(Simulation&, int);

			Simulation& sim;

			double3 pos;
			double3 n;
			double3 s1;
			double3 s2;
#ifdef SPOTS
			vector<double3> spot; // position of spot in space
			vector<double3> spotNorm; // normalized direction of spot
#endif

			int ID;

			void updateSS();
			void updateTo(double3, double3);
			void move(double3);
			void rotate(double3, double, bool);
			void rotate(double3, double, double3);
#ifdef SPOTS
			void moveSpot(int, double3);
			bool spotIsOnCorrectHemisphere(int) const;
			double spotE() const;
			void updateSpotPos();
			void updateSpotNorm();
#endif
			bool intersects() const;
		};

		vector<Spherocyl> spherocyl;
		double beta;
		double E;

		double3 size;

		double scN;
		double scR;
		double scD;
		double scT;
		double scTr;
		double scTfactor;
		double scOneOverT;
		double scOneOverTrad;
		bool scAngleLimited;
		double scTC;
		double sc2R2; // distance for intersection
		double scC2; // distance cutoff for center of spheres
		bool maxTwoInters;

#ifdef SPOTS
		int scSpots;
		double scK; // spot EM weight factor
		double scD2;
		int countInterSpots(const Spherocyl&, const Spherocyl&, bool, bool);
#endif

		// MC moves
		MCfeedback MCrnd();
		MCfeedback MCmove();
		MCfeedback MCrot();
		MCfeedback MCcluMove();
		MCfeedback MCcluRot();
		MCfeedback MCsnake();
		MCfeedback MCcosSnake();
#ifdef SPOTS
		MCfeedback MCspotMove();
#endif
		bool decide(double);

		// VTF stuff
		string VTFfile;
		string SimFile;
		void saveVTF() const;
		void saveVTFheader() const;
		unsigned long VTFfreq = 1'000;
		unsigned long VTFelapsed = 0;

		// random
		Randomizer r;

		// grid things
		vector3<vector<Spherocyl*>> grid;
		int3 gridCells;
		double gridSize;
		void makeGrid();
		int3 getGridPos(double3) const;
		void updateGrid(Spherocyl*, double3, double3);


		// energy calculation things
		void calcE();
		double calcEcontribution(const Spherocyl&);
		double calcInter(const Spherocyl&, const Spherocyl&) const;
		bool intersectsSomeoneIgnoreGrid(const Spherocyl&) const;
		double calcV(double) const;

		// cluster things
		vector<Spherocyl*> getCluster(Spherocyl&);
		vector<Spherocyl*> getInteractingSpherocyls(Spherocyl*);
		void move(vector<Spherocyl*>&, double3);
		void rotate(vector<Spherocyl*>&, double3, double, double3);
		bool intersects(vector<Spherocyl*>&);
		double calcEcontribution(const vector<Spherocyl*>);

		void findOddities() const;
		void saveSimInfo();

		bool equals(const vector<Spherocyl*>&, const vector<Spherocyl*>&) const;

		// snake things
		vector<double> cosValue;
		vector<double> cosProb;
		vector<double> cosProbBracket;

		bool cosLoaded = false;
		vector<pair<double,double>> getBondCoss(Spherocyl&);
		double getProposalProb(Spherocyl&);
		int countBonds(Spherocyl&);

	public:
		Simulation(SimSettings, string);
		IterationFeedback iterate(unsigned long);
		void setBeta(double);
		void setVTFfreq(unsigned long);
		double getE();

		// PBC things
		double3 getPBC(const double3&) const;
		double3 dispPBC(const double3&, const double3&) const;
		double distPBC(const double3& a, const double3& b) const;
		double dist2PBC(const double3& a, const double3& b) const;
		double dist2PBC(const Spherocyl& a, const Spherocyl& b) const;

		double calcVolFrac() const;
		double calcOrderParameterLocal();
		double calcOrderParameterGlobal();
		double calcBiggestClusterFraction();
		double3 calcClusterSizes(double&);
		void updateCosAngleCount(vector<int>&);
		void loadCosDistr(const char* filename);
		void saveInteraction(int, string);

		// scaffolding network stuff
		vector<vector<bool>> getSphereAdjMatrix();
		vector<vector<int>> getSphereClusters(const vector<vector<bool>>&) const;
		vector<vector<int>> getScafAdjMatrix(const vector<vector<bool>>&, const vector<vector<int>>&) const;
		vector<int> getEdgeSizes(const vector<vector<int>>&) const;



};

Simulation::Simulation(SimSettings s, string file)
{
	VTFfile = file + ".vtf";
	SimFile = file + ".info";

	beta = 0.0;
	scN = s.sc;
	scR = s.scR;
	scD = s.scD;
	scOneOverT = s.scOneOverT;
	scOneOverTrad = scOneOverT * 180.0 / 3.14159;
	scTr = 1.0 / scOneOverTrad;
	scT = scTr * 180.0/3.14159;
	scTfactor = scOneOverTrad * 3.14159;
	scTC = cos(scTr);
	scAngleLimited = (scOneOverT > 1.0/180.0);
	sc2R2 = 4*scR*scR;
	scC2 = (2*scR+scD)*(2*scR+scD);
	maxTwoInters = (4*scR+2*scD < 2.0);
	//cerr << maxTwoInters << endl;
	cerr << scT << endl;

#ifdef SPOTS
	scSpots = s.scSpots;
	scK = s.scK;
	scD2 = scD*scD;
				next.spotNorm.push_back(r.randomSemispherical(next.n));
				next.spot.push_back(usePBC(next.pos + next.n + scR * next.spotNorm[i]));
#endif

	gridSize = 2.0+2*(scR+scD/2);
	size = s.size * gridSize;
	gridCells = s.size;

	// allocate memory n stuff, place randomly

	unsigned long trials = 0;
	unsigned long maxTrials = s.sc * 100'000; //100000 trials per stick

	const double critVolFrac = 0.208;
	bool initRandom = (calcVolFrac() < critVolFrac);

	if (initRandom)
	{
		do
		{
			for (unsigned int i = 0; i < scN; i++)
			{
				Spherocyl next(*this, i);

				// make sure it does not intersect
				do
				{
					next.pos = { size.x*r.random(), size.y*r.random(), size.z*r.random() };
					next.n = r.randomSpherical();
					trials++;

					if (trials >= maxTrials) {break;}
				} while (intersectsSomeoneIgnoreGrid(next));

				if (trials >= maxTrials)
				{
					spherocyl.clear();
					break;
				}

#ifdef SPOTS
				for (int i = 0; i < scSpots; i++)
				{
					next.spotNorm.push_back(r.randomSemispherical(next.n));
					next.spot.push_back(usePBC(next.pos + next.n + scR * next.spotNorm[i]));
				}
				for (int i = 0; i < scSpots; i++)
				{
					next.spotNorm.push_back(r.randomSemispherical(-(next.n)));
					next.spot.push_back(usePBC(next.pos - next.n + scR * next.spotNorm[i+scSpots]));
				}
#endif

				next.updateSS();
				spherocyl.push_back(next);

			}
		} while (s.sc != spherocyl.size());
	}
	else
	{
		// initialize nice and toight
		int inX, inY, inZ;
		double dispXY = 2*scR + 1E-3;
		double dispZ = 2*scR+2.0+1E-3;
		inX = (int)(size.x/(dispXY));
		inY = (int)(size.y/dispXY);
		inZ = (int)(size.z/dispZ);

		int maxN = inX*inY*inZ;
		cout << maxN << endl;

		if (scN > maxN)
		{
			cerr << "that's just too much, man!" << endl;
			throw 1;
		}

		int posX, posY, posZ;
		posX = 0;
		posY = 0;
		posZ = 0;

		for (unsigned int i = 0; i < scN; i++)
		{
			Spherocyl next(*this, i);


			next.pos = { posX*dispXY, posY*dispXY, posZ*dispZ};
			next.n = {0,0,1};
			next.updateSS();
			spherocyl.push_back(next);

			posX++;
			if (posX >= inX) {posY++; posX=0;}
			if (posY >= inY) {posZ++; posY=0;}
		}

	}



	// initialize grid
	makeGrid();

	// calculate energy
	calcE();

	// save VTF
	saveVTFheader();
	saveVTF();

	// save info file
	saveSimInfo();
}

#ifdef SPOTS
int Simulation::countInterSpots(const Spherocyl& a, const Spherocyl& b, bool nap, bool nbp)
{
	int toReturn = 0;

	for (int i = (nap?0:scSpots); i < (nap?scSpots:2*scSpots); i++)
	{
		for (int j = (nbp?0:scSpots) + i; j < (nbp?scSpots:scSpots*2); j++)
		{
			if (dist2PBC( a.spot[i], b.spot[j] ) < scD2)
			{
				toReturn++;
			}
		}
	}

	return toReturn;
}
#endif

Simulation::Spherocyl::Spherocyl(Simulation& s, int id) : sim(s), ID(id)
{

}

double Simulation::calcVolFrac() const
{
	double volSphere = 4.0/3.0*PI*scR*scR*scR;
	double areaCirc = PI*scR*scR;
	double volCyl = 2.0*areaCirc;

	return scN*(volSphere+volCyl)/(size.x*size.y*size.z);
}

double Simulation::calcOrderParameterLocal()
{
	double r2 = 4.0*gridSize*gridSize;

	int pairs = 0;
	double sum = 0;

	for (const auto& s : spherocyl)
	{
		int3 gp = getGridPos(s.pos);

		for (int dx = -2; dx <= 2; dx++)
		{
			for (int dy = -2; dy <=2; dy++)
			{
				for (int dz = -2; dz <= 2; dz++)
				{
					int3 pp = { (gp.x+dx+gridCells.x)%gridCells.x,
						(gp.y+dy+gridCells.y)%gridCells.y,
						(gp.z+dz+gridCells.z)%gridCells.z };

					const auto& nv = grid.at(pp);

					for (const auto& partner : nv)
					{
						if (partner -> ID != s.ID && dist2PBC(partner -> pos, s.pos) < r2 )
						{
							double cosTheta =  partner->n * s.n;
							sum += (3.0*cosTheta*cosTheta-1.0)*0.5;
							pairs++;
						}
					}

				}
			}
		}
	}

	return (pairs > 0 ? sum / pairs : 0.0);
}

double Simulation::calcOrderParameterGlobal()
{
	int pairs = 0;
	double sum = 0;

	for (const auto& s : spherocyl)
	{
		for (const auto& partner : spherocyl)
		{
			if (partner.ID > s.ID)
			{
				double cosTheta =  partner.n * s.n;
				sum += (3.0*cosTheta*cosTheta-1.0)*0.5;
				pairs++;
			}
		}

	}

	return (pairs > 0 ? sum / pairs : 0.0);
}

double Simulation::calcBiggestClusterFraction()
{
	vector<bool> foundAlready(scN, false);
	int biggestCluster = 0;

	for (auto&& s : spherocyl)
	{
		if (!foundAlready[s.ID])
		{
			const auto& newCluster = getCluster(s);

			for (const auto& c : newCluster)
			{
				foundAlready[c->ID] = true;
			}

			if (newCluster.size() > biggestCluster) {biggestCluster = newCluster.size();}
		}
	}

	return (double)biggestCluster/scN;
}

double3 Simulation::calcClusterSizes(double& av) // calculates biggest cluster (x), 2nd biggest (y) and reduced cluster size (z)
{
	vector<int> clusterSizes;
	vector<bool> foundAlready(scN, false);

	for (auto&& s : spherocyl)
	{
		if (!foundAlready[s.ID])
		{
			const auto& newCluster = getCluster(s);

			for (const auto& c : newCluster)
			{
				foundAlready[c->ID] = true;
			}

			clusterSizes.push_back(newCluster.size());
		}
	}

	// sort
	sort(clusterSizes.begin(), clusterSizes.end());

	double3 toReturn;
	toReturn.x = clusterSizes[clusterSizes.size() - 1];
	toReturn.y = clusterSizes[clusterSizes.size() - 2];

	// for reduced cluster size
	vector<int> clusterSizeCount;
	vector<int> clusterSize;

	int currentSize = clusterSizes[clusterSizes.size() - 1];
	int currentCount = 0;
	av = 0.0;
	for (const auto& i : clusterSizes)
	{
		av += i;
		if (currentSize == i)
		{
			currentCount++;
		}
		else
		{
			clusterSizeCount.push_back(currentCount);
			clusterSize.push_back(currentSize);

			currentSize = i;
			currentCount = 1;
		}
	}
	av /= clusterSizes.size();

	clusterSizeCount.push_back(currentCount);
	clusterSize.push_back(currentSize);

	int nomin = 0;
	int denomin = 0;

	for (int j = 0; j < clusterSize.size(); j++)
	{
		int i = clusterSize[j];
		int Ni = clusterSizeCount[j];

		nomin += Ni*i*i;
		denomin += Ni*i;
	}

	nomin -= clusterSize[clusterSize.size()-1] * clusterSize[clusterSize.size() - 1] * clusterSizeCount[clusterSize.size() - 1];

	toReturn.z = ((double)nomin)/denomin;

	return toReturn;
}



void Simulation::saveSimInfo()
{
	ofstream file(SimFile.c_str());

	file << SimFile << endl;
	file << "scR: " << scR << endl;
	file << "scD: " << scD << endl;
	file << "scT: " << scT << endl;
	file << "scOneOverT: " << scOneOverT << endl;
	file << "scAngleLimited: " << scAngleLimited << endl;
	file << "maxTwoInters: " << maxTwoInters << endl;
	file << "spherocyl.size(): " << spherocyl.size() << endl;
	file << "size: " << size << endl;
	file << "gridCells: " << gridCells << endl;
	file << "volume fraction: " << calcVolFrac() << endl;
	file << "use spots: ";
#ifdef SPOTS
	file << "yes" << endl;
	file << "spots: " << scSpots << endl;
	file << "spot EM repulsion (scK): " << scK << endl;
#else
	file << "no" << endl;
#endif

	file.close();
}

bool Simulation::equals(const vector<Spherocyl*>& a, const vector<Spherocyl*>& b) const
{
	vector<bool> inA(scN, false);
	vector<bool> inB(scN, false);

	for (const auto& aa : a)
	{
		inA[ aa -> ID ] = true;
	}

	for (const auto& bb : b)
	{
		inB[ bb -> ID ] = true;
	}

	for (int i = 0; i < scN; i++)
	{
		if (inA[i] != inB[i]) { return false; }
	}

	return true;
}

double Simulation::calcV(double c) const
{
	double angle = ( c >= 1.0 ? 0.0 : ( c <= -1.0 ? fp::pi : acos(c) ) );
	/*if (c < -1 || c > 1)
	  {
	  cerr << "c outside range: " << c << endl;
	  }*/

	double cc = cos(angle*scTfactor);

	return .5 + .5*cc;
}

double Simulation::calcInter(const Spherocyl& a, const Spherocyl& b) const
{
	double toReturn = 0.0;

	const double3& a1 = a.s1;
	const double3& a2 = a.s2;
	const double3& b1 = b.s1;
	const double3& b2 = b.s2;
	const double3& na = a.n;
	const double3& nb = b.n;

#ifdef SPOTS
	if (maxTwoInters)
	{
		// only counts correctly if there can be at most one bond between two spherocyl and one sphere
		if (dist2PBC(a1, b1) < scC2)
		{
			toReturn += countInterSpots(a, b, true, true);
			if (dist2PBC(a2, b2) < scC2)
			{
				toReturn += countInterSpots(a, b, false, false);
			}
		}
		else if (dist2PBC(a1, b2) < scC2)
		{
			toReturn += countInterSpots(a, b, true, false);

			if (dist2PBC(a2, b1) < scC2)
			{
				toReturn += countInterSpots(a, b, false, true);
			}
		}
		else if (dist2PBC(a2, b1) < scC2)
		{
			toReturn += countInterSpots(a, b, false, true);
		}
		else if (dist2PBC(a2, b2) < scC2)
		{
			toReturn += countInterSpots(a, b, false, false);
		}
	}
	else
	{
		if (dist2PBC(a1, b1) < scC2)
		{
			toReturn += countInterSpots(a, b, true, true);
		}
		if (dist2PBC(a1, b2) < scC2)
		{
			toReturn += countInterSpots(a, b, true, false);
		}
		if (dist2PBC(a2, b1) < scC2)
		{
			toReturn += countInterSpots(a, b, false, true);
		}
		if (dist2PBC(a2, b2) < scC2)
		{
			toReturn += countInterSpots(a, b, false, false);
		}
	}
#else
	if (maxTwoInters)
	{
		// only counts correctly if there can be at most one bond between two spherocyl and one sphere
		if (dist2PBC(a1, b1) < scC2)
		{
			double3 a1b1 = normalize(dispPBC(a1,b1));
			double c1 = na*a1b1;
			double c2 = -nb*a1b1;

			if (!scAngleLimited || (c1 > scTC && c2 > scTC)) {toReturn -= calcV(c1) * calcV(c2);}

			if (dist2PBC(a2, b2) < scC2)
			{
				double3 a2b2 = normalize(dispPBC(a2,b2));
				double d1 = -na*a2b2;
				double d2 = nb*a2b2;
				if (!scAngleLimited || (d1 > scTC && d2 > scTC)) {toReturn -= calcV(d1) * calcV(d2);}
			}
		}
		else if (dist2PBC(a1, b2) < scC2)
		{
			double3 a1b2 = normalize(dispPBC(a1,b2));
			double c1 = na*a1b2;
			double c2 = nb*a1b2;
			if (!scAngleLimited || (c1 > scTC && c2 > scTC)) {toReturn -= calcV(c1) * calcV(c2);}

			if (dist2PBC(a2, b1) < scC2)
			{
				double3 a2b1 = normalize(dispPBC(a2,b1));
				double d1 = -na*a2b1;
				double d2 = -nb*a2b1;
				if (!scAngleLimited || (d1 > scTC && d2 > scTC)) {toReturn -= calcV(d1)*calcV(d2);}
			}
		}
		else if (dist2PBC(a2, b1) < scC2)
		{
			double3 a2b1 = normalize(dispPBC(a2,b1));
			double c1 = -na*a2b1;
			double c2 = -nb*a2b1;
			if (!scAngleLimited || (c1 > scTC && c2 > scTC)) {toReturn -= calcV(c1)*calcV(c2);}
		}
		else if (dist2PBC(a2, b2) < scC2)
		{
			double3 a2b2 = normalize(dispPBC(a2,b2));
			double c1 = -na*a2b2;
			double c2 = nb*a2b2;
			if (!scAngleLimited || (c1 > scTC && c2 > scTC)) {toReturn -= calcV(c1)*calcV(c2);}
		}
	}
	else
	{
		if (dist2PBC(a1, b1) < scC2)
		{
			double3 a1b1 = normalize(dispPBC(a1,b1));
			double c1 = na*a1b1;
			double c2 = -nb*a1b1;

			if (!scAngleLimited || (c1 > scTC && c2 > scTC)) {toReturn -= calcV(c1) * calcV(c2);}
		}
		if (dist2PBC(a1, b2) < scC2)
		{
			double3 a1b2 = normalize(dispPBC(a1,b2));
			double c1 = na*a1b2;
			double c2 = nb*a1b2;
			if (!scAngleLimited || (c1 > scTC && c2 > scTC)) {toReturn -= calcV(c1) * calcV(c2);}
		}
		if (dist2PBC(a2, b1) < scC2)
		{
			double3 a2b1 = normalize(dispPBC(a2,b1));
			double c1 = -na*a2b1;
			double c2 = -nb*a2b1;
			if (!scAngleLimited || (c1 > scTC && c2 > scTC)) {toReturn -= calcV(c1)*calcV(c2);}
		}
		if (dist2PBC(a2, b2) < scC2)
		{
			double3 a2b2 = normalize(dispPBC(a2,b2));
			double c1 = -na*a2b2;
			double c2 = nb*a2b2;
			if (!scAngleLimited || (c1 > scTC && c2 > scTC)) {toReturn -= calcV(c1)*calcV(c2);}
		}
	}
#endif

	return toReturn;
}

void Simulation::Spherocyl::updateSS()
{
	s1 = sim.getPBC(pos+n);
	s2 = sim.getPBC(pos-n);
}

void Simulation::Spherocyl::updateTo(double3 np, double3 nn)
{
	double3 oldPos = pos;
	pos = np;
	n = nn;
	sim.updateGrid(this, oldPos, pos);
	updateSS();

#ifdef SPOTS
	updateSpotNorm();
	updateSpotPos();
#endif
}

void Simulation::Spherocyl::move(double3 d)
{
	double3 oldPos = pos;
	pos = sim.getPBC(pos + d);
	sim.updateGrid(this, oldPos, pos);
	updateSS();
#ifdef SPOTS
	updateSpotPos();
#endif
}

void Simulation::Spherocyl::rotate(double3 axis, double angle, bool orig)
{
	double3 oldN = n;
	n = fp::rotate(oldN, axis, angle);
	double3 s1disp = n - oldN;
	double3 oldPos = pos;
	if (orig)
	{
		pos = sim.getPBC(pos - s1disp);
	}
	else
	{
		pos = sim.getPBC(pos + s1disp);
	}
	sim.updateGrid(this, oldPos, pos);
	updateSS();
#ifdef SPOTS
	updateSpotNorm();
	updateSpotPos();
#endif
}

void Simulation::Spherocyl::rotate(double3 axis, double angle, double3 orig)
{
	double3 oldPos = pos;

	pos = sim.getPBC(orig + fp::rotate(sim.dispPBC(orig, pos), axis, angle));
	n = fp::rotate(n, axis, angle);
	updateSS();

	sim.updateGrid(this, oldPos, pos);
}

int3 Simulation::getGridPos(double3 d) const
{
	int3 i;

	i.x = (int)(d.x/gridSize);
	i.y = (int)(d.y/gridSize);
	i.z = (int)(d.z/gridSize);

	return i;
}

bool Simulation::Spherocyl::intersects() const
{
	int3 post;

	for (int dx = -1; dx <= 1; dx++)
	{
		for (int dy = -1; dy <= 1; dy++)
		{
			for (int dz = -1; dz <= 1; dz++)
			{
				post = sim.getGridPos(pos);
				post.x = (post.x+dx + sim.gridCells.x)%sim.gridCells.x;
				post.y = (post.y+dy + sim.gridCells.y)%sim.gridCells.y;
				post.z = (post.z+dz + sim.gridCells.z)%sim.gridCells.z;

				const auto& v = sim.grid.at(post);

				for (const auto& t : v)
				{
					if (t != this)
					{
						if (sim.dist2PBC(*this, *t) < sim.sc2R2) {return true;}
					}
				}
			}
		}
	}

	return false;
}

#ifdef SPOTS
double Simulation::Spherocyl::spotE() const
{
	double toReturn;
	int spots = spotNorm.size() / 2;

	for (int i = 0; i < 2 * spots; i++)
	{
		for (int j = i + 1; j < (i < spots ? spots : 2*spots); j++)
		{
			toReturn += sim.scK * dist(spotNorm[i], spotNorm[k]);
		}
	}

	return toReturn;
}

double Simulation::Spherocyl::spotE(int a) const
{
	double toReturn;
	int spots = spotNorm.size() / 2;

	for (int j =  (a < spots ? 0 : spots); j < (a < spots ? spots : 2*spots); j++)
	{
		toReturn += sim.scK * dist(spotNorm[a], spotNorm[k]);
	}

	return toReturn;
}

void Simulation::Spherocyl::updateSpotPos()
{
	for (int i = 0; i < sim.scSpots; i++)
	{
		spot[i] = getPBC(s1 + sim.scR * spotNorm[i]);
	}
	for (int i = sim.scSpots; i < 2 * sim.scSpots; i++)
	{
		spot[i] = getPBC(s2 + sim.scR * spotNorm[i]);
	}
}

void Simulation::Spherocyl::updateSpotNorm()//TODO: rotate
#endif

	bool Simulation::intersectsSomeoneIgnoreGrid(const Spherocyl& s) const
{
	for (const auto& t : spherocyl)
	{
		if (&t != &s)
		{
			if (dist2PBC(t, s) < sc2R2) {return true;}
		}
	}

	return false;
}

bool Simulation::intersects(vector<Spherocyl*>& v)
{
	int3 post;

	for (const auto& s : v)
	{
		for (int dx = -1; dx <= 1; dx++)
		{
			for (int dy = -1; dy <= 1; dy++)
			{
				for (int dz = -1; dz <= 1; dz++)
				{
					post = getGridPos(s->pos);
					post.x = (post.x+dx + gridCells.x)%gridCells.x;
					post.y = (post.y+dy + gridCells.y)%gridCells.y;
					post.z = (post.z+dz + gridCells.z)%gridCells.z;

					auto&& v2 = grid.at(post);

					for (const auto& t : v2)
					{
						if (s != t)
						{
							if (dist2PBC(*s, *t) < sc2R2) {return true;}
						}
					}
				}
			}
		}
	}

	return false;
}



void Simulation::updateGrid(Spherocyl* s, double3 p0, double3 p1)
{
	int3 oldGridPos = getGridPos(p0);
	int3 newGridPos = getGridPos(p1);

	if (oldGridPos != newGridPos)
	{
		// delete from old grid pos
		auto& v = grid.at(oldGridPos);
		auto i = find(v.begin(), v.end(), s);
#ifdef debug
		if (i != v.end())
		{
#endif
			v.erase(i);
#ifdef debug
		}
		else
		{
			cerr << "oh no! Trying to remove from grid, but it's just not there!" << endl;
		}
#endif
		// add to new grid pos
		grid.at(newGridPos).push_back(s);
	}
}

double Simulation::distPBC(const double3& a, const double3& b) const
{
	return sqrt(dist2PBC(a,b));
}

double Simulation::dist2PBC(const double3& a, const double3& b) const
{
	//double3 disp = dispPBC(a,b); // TODO: improve?
	double3 up = {fmax(a.x,b.x), fmax(a.y,b.y), fmax(a.z,b.z)};
	double3 lo = {fmin(a.x,b.x), fmin(a.y,b.y), fmin(a.z,b.z)};

	double dx = fmin(up.x - lo.x, lo.x + size.x - up.x);
	double dy = fmin(up.y - lo.y, lo.y + size.y - up.y);
	double dz = fmin(up.z - lo.z, lo.z + size.z - up.z);

	return dx*dx + dy*dy+dz*dz;
}

double3 Simulation::getPBC(const double3& a) const
{
	double dx = a.x - size.x*floor(a.x/size.x);
	double dy = a.y - size.y*floor(a.y/size.y);
	double dz = a.z - size.z*floor(a.z/size.z);

	return double3{dx,dy,dz};
}

double3 Simulation::dispPBC(const double3& a, const double3& b) const
{
	double3 aa = getPBC(a);
	double3 bb = getPBC(b);
	double3 rij = bb - aa;

	double3 res;

	double FX;
	FX = fabs(rij.x);
	if (FX < size.x-FX)
	{
		res.x = rij.x;
	}
	else
	{
		res.x = rij.x - ((rij.x > 0)?size.x:-size.x);
	}

	FX = fabs(rij.y);
	if (FX < size.y-FX)
	{
		res.y = rij.y;
	}
	else
	{
		res.y = rij.y - ((rij.y > 0)?size.y:-size.y);
	}

	FX = fabs(rij.z);
	if (FX < size.z-FX)
	{
		res.z = rij.z;
	}
	else
	{
		res.z = rij.z - ((rij.z > 0)?size.z:-size.z);
	}

	return res;
}

double Simulation::dist2PBC(const Spherocyl& a, const Spherocyl& b) const
{
	double3 w1 = a.n;
	double3 w2 = b.n;

	double3 min_rij = dispPBC(a.pos, b.pos);

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

IterationFeedback Simulation::iterate(unsigned long iterations)
{
#ifndef NOFEEDBACK
	IterationFeedback toReturn;
#endif
	for (unsigned long iter = 0; iter < iterations; iter++)
	{
		// pick a random move
		double rand = r.random();

		/*if (rand < 0.5)
		  {
		  MCrnd();
		  }
		  else
		  {
		  MCsnake();
		  }*/

		if (rand < 0.2)
		{
#ifndef NOFEEDBACK
			toReturn.move += MCmove().dist;
#else
			MCmove();
#endif
		}
		else if (rand < 0.4)
		{
			MCrnd();
		}
		else if (rand < 0.6)
		{
			if (cosLoaded) {MCcosSnake();} else {MCsnake();}
		}
		else if (rand < 0.98)
		{
#ifndef NOFEEDBACK
			toReturn.rot += MCrot().dist;
#else
			MCrot();
#endif
		}
		else if (rand < 0.99)
		{
#ifndef NOFEEDBACK
			toReturn.cluMove += MCcluMove().dist;
#else
			MCcluMove();
#endif
		}
		else
		{
#ifndef NOFEEDBACK
			toReturn.cluRot += MCcluRot().dist;
#else
			MCcluRot();
#endif
		}

		// save vtf
		VTFelapsed++;
		if (VTFelapsed >= VTFfreq)
		{
			VTFelapsed = 0;
			saveVTF();
		}
	}

	//findOddities();

#ifndef NOFEEDBACK
	return toReturn;
#else
	return;
#endif
}

bool Simulation::decide(double deltaE)
{
	if (deltaE <= 0)
	{
		return true;
	}

	return (r.random() < exp(-beta*deltaE));
}

MCfeedback Simulation::MCrnd()
{
	// pick a spherocylinder
	auto& s = spherocyl[r.random(spherocyl.size())];

	// propose transition
	double3 newPos = {r.random()*size.x, r.random()*size.y, r.random()*size.z};
	double3 newN = r.randomSpherical();
	double3 oldPos = s.pos;
	double3 oldN = s.n;

	// check if splace available
	s.updateTo(newPos, newN);

	if (s.intersects())
	{
		s.updateTo(oldPos, oldN);
#ifndef NOFEEDBACK
		return {0.0};
#else
		return;
#endif
	}

	// calculate energy difference
	double Eafter = calcEcontribution(s);
	s.updateTo(oldPos, oldN);
	double Ebefore = calcEcontribution(s);

	if (decide(Eafter - Ebefore))
	{
		// accept
		s.updateTo(newPos, newN);
		E+=Eafter-Ebefore;
#ifndef NOFEEDBACK
		return {0.0};
#else
		return;
#endif
	}
	else
	{
		// reject: nothing to do
#ifndef NOFEEDBACK
		return {0.0};
#else
		return;
#endif
	}
}

MCfeedback Simulation::MCmove()
{
	// pick a spherocylinder
	auto& s = spherocyl[r.random(spherocyl.size())];

	// propose transition
	double3 trans = MCmoveDist*r.randomSpherical();

	// check if splace available
	s.move(trans);

	if (s.intersects())
	{
		s.move(-trans);
#ifndef NOFEEDBACK
		return {0.0};
#else
		return;
#endif
	}

	// calculate energy difference
	double Eafter = calcEcontribution(s);
	s.move(-trans);
	double Ebefore = calcEcontribution(s);
#ifndef NOFEEDBACK
	double3 oldPos = s.pos;
#endif

	if (decide(Eafter - Ebefore))
	{
		// accept
		s.move(trans);
		E+=Eafter-Ebefore;
#ifndef NOFEEDBACK
		return {norm(dispPBC(oldPos, s.pos))};
#else
		return;
#endif
	}
	else
	{
		// reject: nothing to do
#ifndef NOFEEDBACK
		return {0.0};
#else
		return;
#endif
	}
}

MCfeedback Simulation::MCrot()
{
	// pick a spherocylinder
	auto& s = spherocyl[r.random(spherocyl.size())];
	double3 oldPos = s.pos;
	double3 oldN = s.n;

	// propose axis
	double3 axis = r.randomSpherical();

	// propose angle and position (save one random call)
	double angle = MCrotAngle*(2*r.random() - 1.0);
	bool posSphere = angle > 0;
	angle = fabs(angle);

	// check if place available
	s.rotate(axis, angle, posSphere);
	double3 newPos = s.pos;
	double3 newN = s.n;

	if (s.intersects())
	{
		//s.rotate(axis, -angle, posSphere);
		s.updateTo(oldPos, oldN);
#ifndef NOFEEDBACK
		return {0.0};
#else
		return;
#endif
	}

	// calculate energy difference
	double Eafter = calcEcontribution(s);
	s.updateTo(oldPos, oldN);
	//s.rotate(axis, -angle, posSphere);
	double Ebefore = calcEcontribution(s);
#ifndef NOFEEDBACK
	double3 oldPos = s.pos;
#endif

	if (decide(Eafter - Ebefore))
	{
		// accept
		s.updateTo(newPos, newN);
		//s.rotate(axis, angle, posSphere);
		E+=Eafter-Ebefore;
#ifndef NOFEEDBACK
		return {norm(dispPBC(oldPos, s.pos))};
#else
		return;
#endif
	}
	else
	{
		// reject: nothing to do
#ifndef NOFEEDBACK
		return {0.0};
#else
		return;
#endif
	}
}

MCfeedback Simulation::MCcluMove()
{
	// pick a spherocylinder
	auto& s = spherocyl[r.random(spherocyl.size())];
	auto cluster = getCluster(s);

	// propose transition
	double3 trans = MCcluMoveDist*r.randomSpherical();

	// check if splace available
#ifndef NOFEEDBACK
	vector<double3> oldPos;
	for (const auto& sc : cluster)
	{
		oldPos.push_back(sc -> pos);
	}
#endif
	move(cluster, trans);

	if (intersects(cluster) || getCluster(s).size() != cluster.size())
	{
		move(cluster, -trans);
#ifndef NOFEEDBACK
		return {0.0};
#else
		return;
#endif
	}
#ifndef NOFEEDBACK
	double toReturn = 0.0;
	for (int i = 0; i < cluster.size(); i++)
	{
		toReturn += norm(dispPBC(oldPos[i], cluster[i] -> pos));
	}
	return {toReturn};
#else
	return;
#endif
}
MCfeedback Simulation::MCcluRot()
{
	// pick a spherocylinder
	auto& s = spherocyl[r.random(spherocyl.size())];
	auto cluster = getCluster(s);

	// save old pos and rot
	vector<double3> oldPos;
	vector<double3> oldN;

	oldPos.resize(cluster.size());
	oldN.resize(cluster.size());
	for (int i = 0; i < cluster.size(); i++)
	{
		oldPos[i] = cluster[i] -> pos;
		oldN[i] = cluster[i] -> n;
	}

	// propose axis
	double3 axis = r.randomSpherical();

	// propose angle
	double angle = MCcluRotAngle*r.random();
	double3 orig = s.pos;

	// check if place available
	rotate(cluster, axis, angle, orig);

	if (intersects(cluster) || !equals(getCluster(s), cluster) )
	{
		// fix
		for (int i = 0; i < cluster.size(); i++)
		{
			cluster[i] -> updateTo(oldPos[i], oldN[i]); // to avoid rotation rounding errors
		}
#ifndef NOFEEDBACK
		return {0.0};
#else
		return;
#endif
	}
	else
	{
		// calc E difference
		double newE = calcEcontribution(cluster);
		// fix
		for (int i = 0; i < cluster.size(); i++)
		{
			cluster[i] -> updateTo(oldPos[i], oldN[i]); // to avoid rotation rounding errors
		}
		double oldE = calcEcontribution(cluster);

		if (newE < oldE || r.random() < exp(-beta*(newE-oldE)))
		{
			// accept
			rotate(cluster, axis, angle, orig);
			E += newE - oldE;

#ifndef NOFEEDBACK
			double toReturn = 0.0;
			for (int i = 0; i < cluster.size(); i++)
			{
				toReturn += norm(dispPBC(oldPos[i], cluster[i] -> pos));
			}
			return {toReturn};
#else
			return;
#endif

		}
		else
		{
			// reject
			//
#ifndef NOFEEDBACK
			double toReturn = 0.0;
			for (int i = 0; i < cluster.size(); i++)
			{
				toReturn += norm(dispPBC(oldPos[i], cluster[i] -> pos));
			}
			return {toReturn};
#else
			return;
#endif

		}
	}

}

MCfeedback Simulation::MCsnake()
{
	// pick rod randomly
	auto& toMove = spherocyl[r.random(scN)];

	// check if it is connected to something
	double oldE = calcEcontribution(toMove);
	bool isBonded = (oldE < 0.0);

	if (!isBonded) {
		//cerr << "nb" << endl;

		return;
	}

	// find partner
	Spherocyl* partner;

	do
	{
		partner = &(spherocyl[r.random(scN)]);
	} while (partner == &toMove);
	//cout << "found partner" << endl;

	// partner is now someone else

	// pick a side, any side
	bool posDirection = (r.random() < 0.5);

	double distance = scR+scR+r.random()*scD;
	double3 disp1;
	double3 normDisp;

	do
	{
		normDisp = r.randomSpherical();
		disp1=distance*normDisp;
	} while (scAngleLimited && normDisp*(posDirection?partner->n:-partner->n) < scTC);
	//cout << "found dir1" << endl;

	double3 newN;
	do
	{
		newN = r.randomSpherical();
	} while (scAngleLimited && normDisp*newN < scTC);
	//cout << "found dir2" << endl;

	double3 newPos = getPBC((posDirection ? partner->s1 : partner->s2) + disp1 + newN);

	// temp. move to new pos
	double3 oldPos = toMove.pos;
	double3 oldN = toMove.n;
	int oldBonds = countBonds(toMove);
	toMove.updateTo(newPos, newN);
	int newBonds = countBonds(toMove);
	//if (oldBonds != newBonds) {cout << oldBonds << " " << newBonds << endl;}

	//cout << "moved" << endl;

	if (toMove.intersects())
	{
		//reject
		toMove.updateTo(oldPos, oldN);
		//cerr << "inter" << endl;
		return;
	}

	double newE = calcEcontribution(toMove);
	if (newE >= 0.0)
	{
		//cerr << newE << endl;
	}
	//cout << "newE" << endl;

	double acc = exp(-beta*(newE-oldE))*(double)oldBonds/(double)newBonds;

	if (acc >= 1.0 || r.random() < acc)
	{
		//accept
		E += newE - oldE;
	}
	else
	{
		toMove.updateTo(oldPos, oldN);
		//cerr << "energ" << endl;
		return;
	}
}

MCfeedback Simulation::MCcosSnake()
{

	//cout << "csnake " << VTFelapsed << endl;
	// pick rod randomly
	auto& toMove = spherocyl[r.random(scN)];

	// check if it is connected to something
	double oldE = calcEcontribution(toMove);
	bool isBonded = (oldE < 0.0);

	if (!isBonded) {
		//cerr << "nb" << endl;

		return;
	}

	// find partner
	Spherocyl* partner;

	do
	{
		partner = &(spherocyl[r.random(scN)]);
	} while (partner == &toMove);
	//cout << "found partner" << endl;

	// partner is now someone else

	// pick a side, any side
	bool posDirection = (r.random() < 0.5);

	double distance = scR+scR+r.random()*scD;

	double3 disp1;
	int slotOfAngle1 = r.randomSlot(cosProb);
	double cosOfAngle1 = r.random(cosValue, slotOfAngle1);
	double3 normDisp = r.getWithFixedCos( (posDirection?partner->n:-partner->n), cosOfAngle1);
	disp1 = distance * normDisp;
	//cout << "found dir1" << endl;

	double3 newN;
	int slotOfAngle2 = r.randomSlot(cosProb);
	double cosOfAngle2 = r.random(cosValue, slotOfAngle2);
	newN = r.getWithFixedCos( normDisp, cosOfAngle2);
	//cout << "found dir2" << endl;

	double3 newPos = getPBC((posDirection ? partner->s1 : partner->s2) + disp1 + newN);

	// temp. move to new pos
	double3 oldPos = toMove.pos;
	double3 oldN = toMove.n;

	double reverseProb = getProposalProb(toMove);
	toMove.updateTo(newPos, newN);
	double forwardProb = getProposalProb(toMove);
	//cout << "moved" << endl;

	if (toMove.intersects())
	{
		//reject
		toMove.updateTo(oldPos, oldN);
		//cerr << "inter" << endl;
		return;
	}

	double newE = calcEcontribution(toMove);
	/*if (newE >= 0.0)
	  {
	  cerr << "newE>=0: " << newE << endl;
	  }*/
	//cout << "newE" << endl;

	double accRate = exp(-beta*(newE-oldE)) * (reverseProb/forwardProb);
	//cout << reverseProb << " " << forwardProb << " " << accRate << endl;
	if (r.random() < accRate)
	{
		//accept
		E += newE - oldE;
	}
	else
	{
		toMove.updateTo(oldPos, oldN);
		//cerr << "energ" << endl;
		return;
	}
}

void Simulation::makeGrid()
{
	// allocate grid memory
	grid.resize(gridCells);

	// place spherocylinders on grid
	for (auto&& s : spherocyl)
	{
		int3 gridPos = getGridPos(s.pos);

		grid.at(gridPos).push_back(&s);
	}
}

void Simulation::setBeta(double b)
{
	beta = b;
}

void Simulation::setVTFfreq(unsigned long f)
{
	VTFfreq = f;
	VTFelapsed = 0;
}

double Simulation::getE()
{
	return E;
}

void Simulation::saveVTFheader() const
{
	ofstream file(VTFfile.c_str());

	file << setiosflags(ios::fixed) << setprecision(3);

	// write atoms
	unsigned long atomId = 0;
	for (const auto& s : spherocyl)
	{
		file << "atom " << atomId << " radius " << scR << " name E" << endl;
		atomId++;
		file << "atom " << atomId << " radius " << scR << " name E" << endl;
		atomId++;
	}

	// write bonds
	atomId = 0;
	for (const auto& s : spherocyl)
	{
		file << "bond " << atomId << ":" << atomId+1 << endl;
		atomId+=2;
	}

	file << "pbc " << size.x << " " << size.y << " " << size.z << endl;
	file.close();
}

void Simulation::saveVTF() const
{
	ofstream file(VTFfile.c_str(), ios::app);

	file << setiosflags(ios::fixed) << setprecision(3);

	file << "timestep" << endl;
	file << "pbc " << size.x << " " << size.y << " " << size.z << endl;

	// write atom positions
	unsigned long atomId = 0;
	for (const auto& s : spherocyl)
	{
		double3 s1 = s.pos + s.n;
		double3 s2 = s.pos - s.n;
		file << s1.x << " " << s1.y << " " << s1.z << endl;
		file << s2.x << " " << s2.y << " " << s2.z << endl;
		atomId++;
	}

	file.close();
}

void Simulation::calcE()
{
	E = 0.0;

	int3 pos;

	for (const auto& s0 : spherocyl)
	{
		E += calcEcontribution(s0);
	}

	E /= 2; // remove double counting

#ifdef SPOTS
	// add spot energy
	for (const auto& s0 : spherocyl)
	{
		E += s0.spotE();
	}
#endif
}

double Simulation::calcEcontribution(const Spherocyl& s)
{
	double toReturn = 0.0;
	int3 pos;

	for (int dx = -1; dx <= 1; dx++)
	{
		for (int dy = -1; dy <= 1; dy++)
		{
			for (int dz = -1; dz <= 1; dz++)
			{
				pos = getGridPos(s.pos);
				pos.x = (pos.x + dx + gridCells.x)%gridCells.x;
				pos.y = (pos.y + dy + gridCells.y)%gridCells.y;
				pos.z = (pos.z + dz + gridCells.z)%gridCells.z;

				auto& v = grid.at(pos);

				for (const auto& t : v)
				{
					if (t != &s)
					{

						toReturn += calcInter(s, *t);
					}
				}
			}
		}
	}

	return toReturn;
}

double Simulation::calcEcontribution(const vector<Spherocyl*> v)
{
	// create bool of who is in the cluster
	vector<bool> isInCluster(scN, false);
	for (const auto& a : v) {isInCluster[a->ID]=true;}

	double toReturn = 0.0;
	int3 pos;

	for (const auto& a : v)
	{
		for (int dx = -1; dx <= 1; dx++)
		{
			for (int dy = -1; dy <= 1; dy++)
			{
				for (int dz = -1; dz <= 1; dz++)
				{
					pos = getGridPos(a -> pos);
					pos.x = (pos.x + dx + gridCells.x)%gridCells.x;
					pos.y = (pos.y + dy + gridCells.y)%gridCells.y;
					pos.z = (pos.z + dz + gridCells.z)%gridCells.z;

					auto& v2 = grid.at(pos);

					for (const auto& t : v2)
					{
						if (a != t && (!isInCluster[ t->ID ] || t->ID > a -> ID) ) // avoid double-counting of inter-cluster-energies
						{

							toReturn += calcInter(*a, *t);
						}
					}
				}
			}
		}
	}
	return toReturn;
}

vector<Simulation::Spherocyl*> Simulation::getCluster(Spherocyl& s)
{
	vector<Spherocyl*> toReturn;

	vector<Spherocyl*> foundLastRound;
	vector<Spherocyl*> foundThisRound;
	vector<bool> isPartOfCluster(spherocyl.size(), false);

	isPartOfCluster[s.ID] = true;

	foundLastRound.push_back(&s);
	toReturn.push_back(&s);

	bool foundNewSpherocyl = true;

	while (foundNewSpherocyl)
	{
		foundNewSpherocyl = false;

		for (auto&& b : foundLastRound)
		{
			auto closeBySpherocyls = getInteractingSpherocyls(b);

			for (auto&& c : closeBySpherocyls)
			{
				int newSpID = c -> ID;

				if (!isPartOfCluster[newSpID])
				{
					isPartOfCluster[newSpID] = true;

					foundThisRound.push_back(&spherocyl[newSpID]);

					foundNewSpherocyl = true;
				}
			}
		}

		for (auto& p : foundThisRound)
		{
			toReturn.push_back(p);
		}

		foundLastRound = foundThisRound;
		foundThisRound.clear();
	}

	return toReturn;
}

vector<Simulation::Spherocyl*> Simulation::getInteractingSpherocyls(Spherocyl* s)
{
	vector<Spherocyl*> toReturn;
	int3 pos;

	for (int dx = -1; dx <= 1; dx++)
	{
		for (int dy = -1; dy <= 1; dy++)
		{
			for (int dz = -1; dz <= 1; dz++)
			{
				pos = getGridPos(s->pos);
				pos.x = (pos.x + dx + gridCells.x)%gridCells.x;
				pos.y = (pos.y + dy + gridCells.y)%gridCells.y;
				pos.z = (pos.z + dz + gridCells.z)%gridCells.z;

				auto& v = grid.at(pos);

				for (const auto& t : v)
				{
					if (t != s)
					{
						if (calcInter(*s, *t) < 0) {toReturn.push_back(t);}
					}
				}
			}
		}
	}

	return toReturn;

}

void Simulation::move(vector<Spherocyl*>& c, double3 d)
{
	for (auto&& s : c)
	{
		s -> move(d);
	}
}

void Simulation::rotate(vector<Spherocyl*>& c, double3 axis, double angle, double3 orig)
{
	for (auto&& s : c)
	{
		s -> rotate(axis, angle, orig);
	}
}

void Simulation::findOddities() const
{
	// single sphero stuff
	const double eps = 1e-6;

	for (const auto& s : spherocyl)
	{
		// pos
		if (dist2PBC(s.pos, getPBC(s.pos)) > eps) {cerr << "pos outside" << endl;}

		// length
		if (fabs(norm2(s.n)-1.0) > eps) {cerr << "n wrong length" << endl;}
		if (fabs(dist2PBC(s.s1, s.s2) - 4.0) > eps) { cerr << "s1 s2 dist wrong: " << dist2PBC(s.s1, s.s2) << " " << s.s1 << " " << s.s2 << endl;}

		// see if its in correct grid pos

	}



}

inline int getPosInArray(double p, int N, double amin, double amax)
{
	return min( N-1, (int)( N*(p-amin)/(amax-amin) ) );
}

void Simulation::updateCosAngleCount(vector<int>& v)
{
	double minC = scTC;
	double maxC = 1.0;
	int N = v.size();

	for (auto&& s1 : spherocyl)
	{
		const auto& inters = getInteractingSpherocyls(&s1);

		for (auto&& s2 : inters)
		{
			const double3& a1 = s1.s1;
			const double3& a2 = s1.s2;
			const double3& b1 = s2->s1;
			const double3& b2 = s2->s2;
			const double3& na = s1.n;
			const double3& nb = s2->n;

			if (maxTwoInters)
			{
				// only counts correctly if there can be at most one bond between two spherocyl and one sphere
				if (dist2PBC(a1, b1) < scC2)
				{
					double3 a1b1 = normalize(dispPBC(a1,b1));
					double c1 = na*a1b1;
					double c2 = -nb*a1b1;

					if (!scAngleLimited || (c1 > scTC && c2 > scTC)) {v[getPosInArray(c1, N, minC, maxC)]++; v[getPosInArray(c2, N, minC, maxC)]++;}

					if (dist2PBC(a2, b2) < scC2)
					{
						double3 a2b2 = normalize(dispPBC(a2,b2));
						double d1 = -na*a2b2;
						double d2 = nb*a2b2;
						if (!scAngleLimited || (d1 > scTC && d2 > scTC)) {v[getPosInArray(d1, N, minC, maxC)]++; v[getPosInArray(d2, N, minC, maxC)]++;}
					}
				}
				else if (dist2PBC(a1, b2) < scC2)
				{
					double3 a1b2 = normalize(dispPBC(a1,b2));
					double c1 = na*a1b2;
					double c2 = nb*a1b2;
					if (!scAngleLimited || (c1 > scTC && c2 > scTC)) {v[getPosInArray(c1, N, minC, maxC)]++; v[getPosInArray(c2, N, minC, maxC)]++;}

					if (dist2PBC(a2, b1) < scC2)
					{
						double3 a2b1 = normalize(dispPBC(a2,b1));
						double d1 = -na*a2b1;
						double d2 = -nb*a2b1;
						if (!scAngleLimited || (d1 > scTC && d2 > scTC)) {v[getPosInArray(d1, N, minC, maxC)]++; v[getPosInArray(d2, N, minC, maxC)]++;}
					}
				}
				else if (dist2PBC(a2, b1) < scC2)
				{
					double3 a2b1 = normalize(dispPBC(a2,b1));
					double c1 = -na*a2b1;
					double c2 = -nb*a2b1;
					if (!scAngleLimited || (c1 > scTC && c2 > scTC)) {v[getPosInArray(c1, N, minC, maxC)]++; v[getPosInArray(c2, N, minC, maxC)]++;}
				}
				else if (dist2PBC(a2, b2) < scC2)
				{
					double3 a2b2 = normalize(dispPBC(a2,b2));
					double c1 = -na*a2b2;
					double c2 = nb*a2b2;
					if (!scAngleLimited || (c1 > scTC && c2 > scTC)) {v[getPosInArray(c1, N, minC, maxC)]++; v[getPosInArray(c2, N, minC, maxC)]++;}
				}
			}
			else
			{
				if (dist2PBC(a1, b1) < scC2)
				{
					double3 a1b1 = normalize(dispPBC(a1,b1));
					double c1 = na*a1b1;
					double c2 = -nb*a1b1;

					if (!scAngleLimited || (c1 > scTC && c2 > scTC)) {v[getPosInArray(c1, N, minC, maxC)]++; v[getPosInArray(c2, N, minC, maxC)]++;}
				}
				if (dist2PBC(a1, b2) < scC2)
				{
					double3 a1b2 = normalize(dispPBC(a1,b2));
					double c1 = na*a1b2;
					double c2 = nb*a1b2;
					if (!scAngleLimited || (c1 > scTC && c2 > scTC)) {v[getPosInArray(c1, N, minC, maxC)]++; v[getPosInArray(c2, N, minC, maxC)]++;}
				}
				if (dist2PBC(a2, b1) < scC2)
				{
					double3 a2b1 = normalize(dispPBC(a2,b1));
					double c1 = -na*a2b1;
					double c2 = -nb*a2b1;
					if (!scAngleLimited || (c1 > scTC && c2 > scTC)) {v[getPosInArray(c1, N, minC, maxC)]++; v[getPosInArray(c2, N, minC, maxC)]++;}
				}
				if (dist2PBC(a2, b2) < scC2)
				{
					double3 a2b2 = normalize(dispPBC(a2,b2));
					double c1 = -na*a2b2;
					double c2 = nb*a2b2;
					if (!scAngleLimited || (c1 > scTC && c2 > scTC)) {v[getPosInArray(c1, N, minC, maxC)]++; v[getPosInArray(c2, N, minC, maxC)]++;}
				}
			}

		}
	}
}


void Simulation::loadCosDistr(const char* filename)
{
	vector<int> count;
	ifstream file(filename);
	int buffer;
	int allCount = 0;
	while (file >> buffer)
	{
		count.push_back(buffer);
		allCount += buffer;
	}
	file.close();

	cosValue.clear();
	cosProb.clear();
	cosProbBracket.clear();
	cosValue.push_back(scTC);
	cosProb.push_back(0.0);
	cosProbBracket.push_back(0.0);

	double binSize = (1.0-scTC)/count.size();

	for (int i = 0; i < count.size(); i++)
	{
		cosProb.push_back( cosProb[i] + (double)(count[i])/allCount );
		cosProbBracket.push_back( (double)(count[i])/allCount );
		cosValue.push_back(cosValue[0] + (i+1)*binSize);
	}
	cosLoaded = true;

	/*for (int i = 0; i < cosValue.size(); i++)
	  {
	  cout << cosValue[i] << " " << cosProb[i] << endl;
	  }*/
}

int Simulation::countBonds(Spherocyl& s)
{
	return getBondCoss(s).size();
}

vector<pair<double,double>> Simulation::getBondCoss(Spherocyl& s1)
{
	vector<pair<double,double>> toReturn;

	const auto& inters = getInteractingSpherocyls(&s1);
	for (auto&& s2 : inters)
	{
		const double3& a1 = s1.s1;
		const double3& a2 = s1.s2;
		const double3& b1 = s2->s1;
		const double3& b2 = s2->s2;
		const double3& na = s1.n;
		const double3& nb = s2->n;

		if (maxTwoInters)
		{
			// only counts correctly if there can be at most one bond between two spherocyl and one sphere
			if (dist2PBC(a1, b1) < scC2)
			{
				double3 a1b1 = normalize(dispPBC(a1,b1));
				double c1 = na*a1b1;
				double c2 = -nb*a1b1;

				if (!scAngleLimited || (c1 > scTC && c2 > scTC))
				{
					toReturn.push_back({c1,c2});
				}

				if (dist2PBC(a2, b2) < scC2)
				{
					double3 a2b2 = normalize(dispPBC(a2,b2));
					double d1 = -na*a2b2;
					double d2 = nb*a2b2;
					if (!scAngleLimited || (d1 > scTC && d2 > scTC))
					{
						toReturn.push_back({d1,d2});
					}
				}
			}
			else if (dist2PBC(a1, b2) < scC2)
			{
				double3 a1b2 = normalize(dispPBC(a1,b2));
				double c1 = na*a1b2;
				double c2 = nb*a1b2;
				if (!scAngleLimited || (c1 > scTC && c2 > scTC))
				{
					toReturn.push_back({c1,c2});
				}

				if (dist2PBC(a2, b1) < scC2)
				{
					double3 a2b1 = normalize(dispPBC(a2,b1));
					double d1 = -na*a2b1;
					double d2 = -nb*a2b1;
					if (!scAngleLimited || (d1 > scTC && d2 > scTC))
					{
						toReturn.push_back({d1,d2});
					}
				}
			}
			else if (dist2PBC(a2, b1) < scC2)
			{
				double3 a2b1 = normalize(dispPBC(a2,b1));
				double c1 = -na*a2b1;
				double c2 = -nb*a2b1;
				if (!scAngleLimited || (c1 > scTC && c2 > scTC))
				{
					toReturn.push_back({c1,c2});
				}
			}
			else if (dist2PBC(a2, b2) < scC2)
			{
				double3 a2b2 = normalize(dispPBC(a2,b2));
				double c1 = -na*a2b2;
				double c2 = nb*a2b2;
				if (!scAngleLimited || (c1 > scTC && c2 > scTC))
				{
					toReturn.push_back({c1,c2});
				}
			}
		}
		else
		{
			if (dist2PBC(a1, b1) < scC2)
			{
				double3 a1b1 = normalize(dispPBC(a1,b1));
				double c1 = na*a1b1;
				double c2 = -nb*a1b1;

				if (!scAngleLimited || (c1 > scTC && c2 > scTC))
				{
					toReturn.push_back({c1,c2});
				}
			}
			if (dist2PBC(a1, b2) < scC2)
			{
				double3 a1b2 = normalize(dispPBC(a1,b2));
				double c1 = na*a1b2;
				double c2 = nb*a1b2;
				if (!scAngleLimited || (c1 > scTC && c2 > scTC))
				{
					toReturn.push_back({c1,c2});
				}
			}
			if (dist2PBC(a2, b1) < scC2)
			{
				double3 a2b1 = normalize(dispPBC(a2,b1));
				double c1 = -na*a2b1;
				double c2 = -nb*a2b1;
				if (!scAngleLimited || (c1 > scTC && c2 > scTC))
				{
					toReturn.push_back({c1,c2});
				}
			}
			if (dist2PBC(a2, b2) < scC2)
			{
				double3 a2b2 = normalize(dispPBC(a2,b2));
				double c1 = -na*a2b2;
				double c2 = nb*a2b2;
				if (!scAngleLimited || (c1 > scTC && c2 > scTC))
				{
					toReturn.push_back({c1,c2});
				}
			}
		}
	}
	return toReturn;
}

double Simulation::getProposalProb(Spherocyl& s)
{
	double toReturn = 0.0;

	//cout << "getting bonds..." << endl;
	auto vp = getBondCoss(s);
	//cout << "got bonds" << endl;

	double c1, c2;

	//cout << "calcing..." << endl;
	for (const auto& p : vp)
	{
		c1 = p.first;
		c2 = p.second;

		double pc1 = cosProbBracket[ r.getSlot(c1, cosValue) ];
		double pc2 = cosProbBracket[ r.getSlot(c2, cosValue) ];

		toReturn += pc1*pc2;
	}
	//cout << "calced" << endl;

	return toReturn;
}

void Simulation::saveInteraction(int type, string filename)
{
	ofstream file(filename.c_str());

	VTFfile = filename+".vtf";
	saveVTFheader();

	// set up system
	spherocyl[0].updateTo({0,0,0},{0,0,1});
	spherocyl[1].updateTo({0,0,2},{0,0,1});

	const int points = 1000;
	for (int i = 0; i < points; i++)
	{
		double3 disp = {0,0,scR+scR+scD/2};
		double3 axis = {1,0,0};
		double angle1 = i * 180.0/points;
		double angle2;

		disp = fp::rotate(disp, axis, i*180.0/points);
		double3 n = {0,0,1};
		if (type==0)
		{
			angle2 = 0;
		}
		else if (type == 1)
		{
			angle2 = 2*angle1;
		}
		else if (type == 2)
		{
			angle2 = angle1;
		}

		n = fp::rotate(n, axis, angle2);

		spherocyl[1].updateTo(spherocyl[0].s1 + disp + n, n);

		if (!spherocyl[0].intersects())
		{
			file << angle1 << " " << angle2 << " " << angleBetween(spherocyl[0].n, spherocyl[1].n) << " " << calcInter(spherocyl[0], spherocyl[1]) << endl;
			VTFfile = filename+".vtf";
			saveVTF();
		}
	}
	file.close();
}

vector<vector<bool>> Simulation::getSphereAdjMatrix()
{
	// allocate
	vector<vector<bool>> toReturn(2*scN);
	for (auto&& i : toReturn)
	{
		i.resize(2*scN);

		for (auto&& j : i) {j = false;}
	}

	// calc matrix
	for (auto&& s1 : spherocyl)
	{
		// diagonal
		toReturn[s1.ID][s1.ID] = true;
		toReturn[s1.ID+1][s1.ID+1] = true;
		// get connections
		const auto& partners = getInteractingSpherocyls(&s1);
		for (auto&& s2 : partners)
		{
			const double3& a1 = s1.s1;
			const double3& a2 = s1.s2;
			const double3& b1 = s2->s1;
			const double3& b2 = s2->s2;
			const double3& na = s1.n;
			const double3& nb = s2->n;

			if (dist2PBC(a1, b1) < scC2)
			{
				double3 a1b1 = normalize(dispPBC(a1,b1));
				double c1 = na*a1b1;
				double c2 = -nb*a1b1;

				if (!scAngleLimited || (c1 > scTC && c2 > scTC))
				{
					toReturn[2*s1.ID][2*s2->ID] = true;
				}
			}
			if (dist2PBC(a1, b2) < scC2)
			{
				double3 a1b2 = normalize(dispPBC(a1,b2));
				double c1 = na*a1b2;
				double c2 = nb*a1b2;
				if (!scAngleLimited || (c1 > scTC && c2 > scTC))
				{
					toReturn[2*s1.ID][2*s2->ID+1] = true;
				}
			}
			if (dist2PBC(a2, b1) < scC2)
			{
				double3 a2b1 = normalize(dispPBC(a2,b1));
				double c1 = -na*a2b1;
				double c2 = -nb*a2b1;
				if (!scAngleLimited || (c1 > scTC && c2 > scTC))
				{
					toReturn[2*s1.ID+1][2*s2->ID] = true;
				}
			}
			if (dist2PBC(a2, b2) < scC2)
			{
				double3 a2b2 = normalize(dispPBC(a2,b2));
				double c1 = -na*a2b2;
				double c2 = nb*a2b2;
				if (!scAngleLimited || (c1 > scTC && c2 > scTC))
				{
					toReturn[2*s1.ID+1][2*s2->ID+1] = true;
				}
			}
		}
	}

	return toReturn;
}

vector<vector<int>> Simulation::getSphereClusters(const vector<vector<bool>>& A) const
{
	vector<vector<int>> toReturn;
	vector<bool> isInACluster(2*scN, false);
	vector<int> currentCluster;

	vector<int> toCheckNextRound;
	vector<int> toCheckThisRound;

	for (int i = 0; i < 2*scN; i++)
	{
		currentCluster.clear();
		toCheckNextRound.clear();
		toCheckThisRound.clear();

		if (!isInACluster[i])
		{
			currentCluster.push_back(i);

			toCheckNextRound.push_back(i);
			isInACluster[i] = true;

			while (toCheckNextRound.size() > 0)
			{
				toCheckThisRound = toCheckNextRound;
				toCheckNextRound.clear();

				for (const auto& j : toCheckThisRound)
				{
					for (int k = 0; k < 2*scN; k++)
					{
						if (A[j][k] && !isInACluster[k])
						{
							isInACluster[k] = true;
							toCheckNextRound.push_back(k);
							currentCluster.push_back(k);
						}
					}
				}
			}

			toReturn.push_back(currentCluster);
		}
	}

	return toReturn;
}

vector<vector<int>> Simulation::getScafAdjMatrix(const vector<vector<bool>>& A, const vector<vector<int>>& c) const
{
	// initialize
	int M = c.size();
	vector<vector<int>> toReturn(M);
	for (auto&& v : toReturn)
	{
		v.resize(M);
		for (auto&& i : v) {i = 0;}
	}

	// lambda to find which cluster has which sphere
	auto cOfs = [&](int sphereId) -> int
	{
		for (int cPos = 0; cPos < M; cPos++)
		{
			const auto& cluster = c[cPos];
			if (find(cluster.begin(), cluster.end(), sphereId) != cluster.end()) {return cPos;}
		}

		cerr << "sphereId " << sphereId << " not found" << endl;
		return -1;
	};

	// fill matrix diagonals
	for (int i = 0; i < M; i++)
	{
		toReturn[i][i] = 1;
	}

	// fill rest of matrix
	for (int s = 0; s < 2*scN; s+=2)
	{
		int c1 = cOfs(s);
		int c2 = cOfs(s+1);

		toReturn[c1][c2] += 1;
		toReturn[c2][c1] += 1;
	}

	return toReturn;
}

vector<int> Simulation::getEdgeSizes(const vector<vector<int>>& A) const
{
	int M = A.size();
	vector<int> toReturn;

	for (int i = 0; i < M; i++)
	{
		for (int j = i+1; j < M; j++)
		{
			if (A[i][j] > 0)
			{
				toReturn.push_back(A[i][j]);
			}
		}
	}

	return toReturn;
}



int main(int args, char** arg)
{
	const long stepsPerIteration = 100'000;
	int N;
	double R;
	double E;
	double T;
	int x,y,z;
	double betaMax;

	if (args <= 11)
	{
		cerr << "Usage: dist [options]...\n" << endl;
		return 1;
	}

	N = atoi(arg[1]);
	R = atof(arg[2]);
	E = atof(arg[3]);
	T = atof(arg[4]);
	x = atoi(arg[5]);
	y = atoi(arg[6]);
	z = atoi(arg[7]);
	betaMax = atof(arg[8]);
	long iterations = atoi(arg[9]);
	int vtfFrames = atoi(arg[10]);
	string fname(arg[11]);

	if (N == -1)
	{
		cout << "saving interaction stuffs" << endl;

		SimSettings settings{2, R, E, T, {x,y,z}};

		Simulation s{settings, (fname).c_str()};

		s.saveInteraction(0, fname+"_parallel");
		s.saveInteraction(1, fname+"_bent");
		s.saveInteraction(2, fname+"_halfbent");

		return 0;
	}

	SimSettings settings{(unsigned int)N, R, E, T, {x,y,z}};

	Simulation s{settings, (fname).c_str()};

	ofstream file((fname + ".data").c_str());
	ofstream histNodeFile((fname + ".hnode").c_str());
	ofstream histEdgeFile((fname + ".hedge").c_str());
	ofstream histCosFile((fname + ".hcos").c_str());

	long vtfFreq = (long)( (double)((long)iterations*(long)stepsPerIteration) / vtfFrames ) ;

	if (args > 12)
	{
		s.loadCosDistr(arg[12]);
	}

	//cerr << "vtfFreq: " << vtfFreq << endl;
	//return 0;

	s.setVTFfreq( vtfFreq ) ;
	for (int it = 1; it <= iterations; it++)
	{
		double beta = betaMax*(double)it/iterations;
		s.setBeta(beta);
		s.iterate(stepsPerIteration);
		double avClusterSize;
		auto CS = s.calcClusterSizes(avClusterSize);

		auto sphereA = s.getSphereAdjMatrix();

		auto sphereClusters = s.getSphereClusters(sphereA);
		auto scafA = s.getScafAdjMatrix(sphereA, sphereClusters);
		auto edgeSizes = s.getEdgeSizes(scafA);
		// calculate average edge numbers
		double avEdge = 0.0;
		for (const auto& e : edgeSizes)
		{
			avEdge += (double)e/edgeSizes.size();
		}

		sample<int> edgeSample(edgeSizes);

		auto hEdge = histogram(edgeSample);

		// save edge histo
		for (int i = 0; i < hEdge.getN(); i++)
		{
			histEdgeFile << it << " " << beta << " " << hEdge[i].pos << " " << hEdge[i].count << endl;
		}


		sample<int> nodeSample;
		for (const auto& i : sphereClusters)
		{
			nodeSample.add(i.size());
		}
		double avNode = 0.0;
		for (const auto& i : sphereClusters)
		{
			avNode += (double)i.size()/sphereClusters.size();
		}

		auto hNode = histogram(nodeSample);
		// save histo
		for (int i = 0; i < hNode.getN(); i++)
		{
			histNodeFile << it << " " << beta << " " << hNode[i].pos << " " << hNode[i].count << endl;
		}


		/*vector<int> his(100,0);
		s.updateCosAngleCount(his);
		histCosFile << it << " " << beta << " ";

		for (const auto& i : his)
		{
			histCosFile << i << " ";
		}
		histCosFile << endl;*/
		file << beta << " " << s.getE() << " " << s.calcOrderParameterLocal() << " " << s.calcOrderParameterGlobal() << " " << CS.x << " " << CS.y << " " << CS.z << " " << (double)CS.x / N << " " << (double) CS.y / N << " " << (double)CS.z / N << " " << avClusterSize << " " << avNode << " " << avEdge << endl;
	}



	file.close();
	histNodeFile.close();
	histEdgeFile.close();
	histCosFile.close();

	return 0;
}
