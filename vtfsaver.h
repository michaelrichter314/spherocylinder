#pragma once

#include "simulation.h"
#include <fstream>

using namespace std;
using namespace fp;

class VTFsaver
{
	public:
		VTFsaver(Simulation&, string);
		void save();

	private:
		string filename;

		void writeHeader();

		Simulation& s;
};

