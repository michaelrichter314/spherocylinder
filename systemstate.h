#pragma once
#include <vector>

// for copying/loading/saving system states
template <typename postype, typename Particletype, typename extrainfo>
struct SystemState
{
	// should contain at least particle positions
	std::vector<Particletype> particle;

	// vector of list of particle IDs, undirected bonds
	// e.g. {{1}, {2,3}, {}, {0}} is bond 0-1, 1-2, 1-3, 3-0
	std::vector<std::vector<int>> bond;

	// system dimension
	postype size;

	// thermodynamical 
	double E;
	double beta;
	double P;

	// additional info, e.g. size and shape of particles
	extrainfo extra;
};
