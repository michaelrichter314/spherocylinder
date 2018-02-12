#pragma once

#include "simulation.h"
#include <fstream>
#include "fpvtfhandler.h"

using namespace std;
using namespace fp;

class VTFloader
{
	public:
		VTFloader(string);
		//SystemState<GridSize, Spherocyl, SpherocylParameter> getState();
		SystemState<GridSize, Spherocyl, SpherocylParameter> getState(int = -1);

		VTFhandler* getHandler();



	private:
		string filename;
		unique_ptr<VTFhandler> handler;

};
