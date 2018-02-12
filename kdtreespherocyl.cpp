#include "kdtreespherocyl.h"

void KdTreeSpherocyl::initialize(const double3& size, Simulation* s)
{
	clear();
	spaceSize = size;

	com.initialize(size, s);
	tip.initialize(size, s);
}

void KdTreeSpherocyl::insert(Spherocyl& s)
{
	const int& id = s.getID();
	com.insert(id, s.getP());
	tip.insert(2*id+0, getPBC(s.getS1()));
	tip.insert(2*id+1, getPBC(s.getS2()));
}

void KdTreeSpherocyl::remove(Spherocyl& s)
{
	const int& id = s.getID();
	com.remove(id, s.getP());
	tip.remove(2*id+0, getPBC(s.getS1()));
	tip.remove(2*id+1, getPBC(s.getS2()));
}

void KdTreeSpherocyl::clear()
{
	com.clear();
	tip.clear();
}

void KdTreeSpherocyl::traverse(TraverseFunction f)
{
	com.traverse(f);
}


bool KdTreeSpherocyl::traverse(TraverseFunction f, const double3& a, const double3& b)
{
	return com.traverse(f, a, b);
}



vector<int> KdTreeSpherocyl::getInBoxCOM(const double3& ll, const double3& ur)
{
	return com.getInBox(ll, ur);
}

vector<int> KdTreeSpherocyl::getInBoxTip(const double3& ll, const double3& ur)
{
	return tip.getInBox(ll, ur);
}

void KdTreeSpherocyl::changeVolume(double3 newSize)
{
	spaceSize = newSize;
	com.changeVolume(newSize);
	tip.changeVolume(newSize);
}

