#include "vtfloader.h"

SystemState<GridSize, Spherocyl, SpherocylParameter> VTFloader::getState(int timestep)
{
	VTFhandler& h = *(handler.get());

	auto t = h.getTimestep(timestep);

	SystemState<GridSize, Spherocyl, SpherocylParameter> toReturn;

	for (int i = 0; i < t.dpos.size(); i+=2)
	{
		double3 p1 = t.dpos[i];
		double3 p2 = t.dpos[i+1];

		Spherocyl sp;
		sp.setS(p1, p2);
		sp.setID(i/2);

		toReturn.particle.push_back(sp);
	}

	// TODO: read additional data
	toReturn.E = 0.0;
	toReturn.beta = 0.0;
	toReturn.P = 0.0;
	toReturn.extra.r = 0.2;
	// attempt to set radius to value from file
	auto vatom = h.getAtoms()[0];
	if (vatom.radius.hasData())
	{
		toReturn.extra.r = vatom.radius.getData();
	}
	toReturn.extra.t = 30;
	toReturn.extra.d = 0.1;

#ifdef grid
	toReturn.size.gridSize = {10,10,10};
	toReturn.size.spaceSize = t.PBC.dsize;
#else
	toReturn.size = t.PBC.dsize;
#endif

	return toReturn;
}

VTFloader::VTFloader(string s) : filename(s)
{
	handler = make_unique<VTFhandler>(s);
};

VTFhandler* VTFloader::getHandler()
{
	return handler.get();
}

