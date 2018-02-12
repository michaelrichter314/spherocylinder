#include "simulation.h"
#include "icmaker.h"
#include "vtfsaver.h"

#include <iostream>

using namespace std;
using namespace fp;

int main(int args, char** arg)
{
	if (args < 6)
	{
		cerr << "usage: " << arg[0] << " <N> <theta> <beta max> <pressure max> <vtf filename>" << endl;
		return 1;
	}

	// read parameters
	int scN = atoi(arg[1]);
	double scT = atof(arg[2]);
	double betaMax = atof(arg[3]);
	double pressureMax = atof(arg[4]);

	// for MC moves and randomizer
	McMover mover;
	
	// make initial condition
	ICmaker icmaker(mover.getRandomizer());

	icmaker.setN(scN);
	icmaker.setSpherocylParameters({0.2, scT, 0.1});

	// set up simulation
	Simulation s(mover.getRandomizer());
	s.setState(icmaker.getState());

	// set up VTF output
	VTFsaver vtf(s, string(arg[5]));

	// save IC
	vtf.save();

	// set up mc mover
	McMove MCSingleDispRot;
	MCSingleDispRot.setFunction(&Simulation::MCSingleGlobalTranslationRotation);
	MCSingleDispRot.setSimulation(&s);
	MCSingleDispRot.enable();
	MCSingleDispRot.setChance(1.0);

	McMove MCSingleRot;
	MCSingleRot.setFunction(&Simulation::MCSingleRotation);
	MCSingleRot.setSimulation(&s);
	MCSingleRot.enable();
	MCSingleRot.setChance(1.0);

	McMove MCSingleDisp;
	MCSingleDisp.setFunction(&Simulation::MCSingleGlobalTranslation);
	MCSingleDisp.setSimulation(&s);
	MCSingleDisp.enable();
	MCSingleDisp.setChance(1.0);

	McMove MCSingleLocalDisp;
	MCSingleLocalDisp.setFunction(&Simulation::MCSingleLocalTranslation);
	MCSingleLocalDisp.setSimulation(&s);
	MCSingleLocalDisp.enable();
	MCSingleLocalDisp.setChance(10.0);

	McMove MCSingleLocalRot;
	MCSingleLocalRot.setFunction(&Simulation::MCSingleLocalRotation);
	MCSingleLocalRot.setSimulation(&s);
	MCSingleLocalRot.enable();
	MCSingleLocalRot.setChance(10.0);

	McMove MCClusterMove;
	MCClusterMove.setFunction(&Simulation::MCClusterMove);
	MCClusterMove.setSimulation(&s);
	MCClusterMove.enable();
	MCClusterMove.setChance(1.0);

#ifdef NPT
	McMove MCVolumeIso;
	MCVolumeIso.setFunction(&Simulation::MCVolumeIso);
	MCVolumeIso.setSimulation(&s);
	MCVolumeIso.enable();
	MCVolumeIso.setChance(0.1);
#endif

	mover.addMove(MCSingleDispRot);
	mover.addMove(MCSingleRot);
	mover.addMove(MCSingleDisp);
	mover.addMove(MCSingleLocalDisp);
	mover.addMove(MCSingleLocalRot);
	mover.addMove(MCClusterMove);
#ifdef NPT
	mover.addMove(MCVolumeIso);
#endif
	const double startBeta = 0.1;
	const double startPressure = 0.1;



	//s.checkSpace();
	int iterations = 1;
	cout << "iterating..." << endl;
	for (int i = 0; i < iterations; i++)
	{
		double beta = startBeta + i*betaMax/iterations;
		double pressure = startPressure + i*pressureMax / iterations;

		s.setBeta(beta);
#ifdef NPT
		s.setPressure(pressure);
#endif
		mover.execute(1'000'000);
		//s.checkSpace();
		cout << i << " "  << beta << " " << pressure << " " << s.getE() << " " << s.getV() << " ";
		for (int i = 0; i < 7; i++)
		{
			cout << mover.getSuccessRate(i) << " ";
		}
		cout << endl;
		vtf.save();
	}
}
