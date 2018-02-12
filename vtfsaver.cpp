#include "vtfsaver.h"

VTFsaver::VTFsaver(Simulation& S, string Filename) : s(S), filename(Filename)
{
	writeHeader();
}

void VTFsaver::save()
{
	ofstream file(filename.c_str(), ofstream::app);

	file << "timestep" << endl;

	SystemState<GridSize, Spherocyl, SpherocylParameter> state = s.getState();

#ifdef grid
	file << "PBC " << state.size.spaceSize.x << " " << state.size.spaceSize.y << " " << state.size.spaceSize.z << endl;
#else
	file << "PBC " << state.size.x << " " << state.size.y << " " << state.size.z << endl;
#endif

	for (auto&& sp : state.particle)
	{
		auto p = sp.getP();
		auto n = sp.getN();

		// calculate s1,s2 manually, so spherocylinders are not cut by PBC
		auto s1 = p+n;
		auto s2 = p-n;

		file << s1.x << " " << s1.y << " " << s1.z << endl;
		file << s2.x << " " << s2.y << " " << s2.z << endl;
	}

	file.close();
}

void VTFsaver::writeHeader()
{
	ofstream file(filename.c_str());

	SystemState<GridSize, Spherocyl, SpherocylParameter> state = s.getState();

	int N = state.particle.size();

	int i = 0;
	for (auto&& sp : state.particle)
	{
		file << "atom " << i << " radius " << state.extra.r << " name S1" << endl;
		i++;
		file << "atom " << i << " radius " << state.extra.r << " name S2" << endl;
		i++;
	}

	i = 0;
	for (auto&& sp : state.particle)
	{
		file << "bond " << i << ":" << i+1<< endl;
		i+=2;
	}

}





