#include "simulation.h"
#include "vtfloader.h"
#include "vtfsaver.h"

#include <iostream>

using namespace std;
using namespace fp;

int main(int args, char** arg)
{
	if (args < 8)
	{
		cerr << "usage: " << arg[0] << " <initial_condition.vtf> <theta> <beta min> <beta max> <pressure min>  <pressure max> <vtf filename>" << endl;
		return 1;
	}

	// read parameters
	string icfile(arg[1]);
	double scT = atof(arg[2]);
	double betaMin = atof(arg[3]);
	double betaMax = atof(arg[4]);
	double pressureMin = atof(arg[5]);
	double pressureMax = atof(arg[6]);

	// for MC moves and randomizer
	McMover mover;
	
	// load initial condition
	VTFloader l(icfile);

	// set up simulation
	Simulation s(mover.getRandomizer());

	auto icstate = l.getState();
	icstate.extra.r = 0.2;
	icstate.extra.d = 0.1;
	icstate.extra.t = scT;

	s.setState(icstate);

	// set up VTF output
	VTFsaver vtf(s, string(arg[7]));

	// save IC
	vtf.save();

	// set up mc mover
	McMove MCSingleDispRot;
	MCSingleDispRot.setFunction(&Simulation::MCSingleGlobalTranslationRotation);
	MCSingleDispRot.setSimulation(&s);
	MCSingleDispRot.enable();
	MCSingleDispRot.setChance(1.0);
	MCSingleDispRot.setImportance(3.0);

	McMove MCSingleRot;
	MCSingleRot.setFunction(&Simulation::MCSingleRotation);
	MCSingleRot.setSimulation(&s);
	MCSingleRot.enable();
	MCSingleRot.setChance(1.0);
	MCSingleRot.setImportance(1.0);

	McMove MCSingleDisp;
	MCSingleDisp.setFunction(&Simulation::MCSingleGlobalTranslation);
	MCSingleDisp.setSimulation(&s);
	MCSingleDisp.enable();
	MCSingleDisp.setChance(1.0);
	MCSingleDisp.setImportance(2.0);

	McMove MCSingleLocalDisp;
	MCSingleLocalDisp.setFunction(&Simulation::MCSingleLocalTranslation);
	MCSingleLocalDisp.setSimulation(&s);
	MCSingleLocalDisp.enable();
	MCSingleLocalDisp.setChance(10.0);
	MCSingleLocalDisp.setImportance(0.3);

	McMove MCSingleLocalRot;
	MCSingleLocalRot.setFunction(&Simulation::MCSingleLocalRotation);
	MCSingleLocalRot.setSimulation(&s);
	MCSingleLocalRot.enable();
	MCSingleLocalRot.setChance(10.0);
	MCSingleLocalRot.setImportance(0.5);

	McMove MCSingleNormal;
	MCSingleNormal.setFunction(&Simulation::MCSingleNormalTranslation);
	MCSingleNormal.setSimulation(&s);
	MCSingleNormal.enable();
	MCSingleNormal.setChance(5.0);
	MCSingleNormal.setImportance(0.1);

	McMove MCClusterMove;
	MCClusterMove.setFunction(&Simulation::MCClusterMove);
	MCClusterMove.setSimulation(&s);
	MCClusterMove.enable();
	MCClusterMove.setChance(1.0);
	MCClusterMove.setImportance(icstate.particle.size());

#ifdef NPT
	McMove MCVolumeIso;
	MCVolumeIso.setFunction(&Simulation::MCVolumeIso);
	MCVolumeIso.setSimulation(&s);
	MCVolumeIso.enable();
	MCVolumeIso.setChance(0.1);
	MCVolumeIso.setImportance(icstate.particle.size());
	//fix chance
	MCVolumeIso.setMinChance(0.1);
	MCVolumeIso.setMaxChance(0.1);

	McMove MCVolumeAniso;
	MCVolumeAniso.setFunction(&Simulation::MCVolumeAniso);
	MCVolumeAniso.setSimulation(&s);
	MCVolumeAniso.enable();
	MCVolumeAniso.setChance(0.1);
	MCVolumeAniso.setImportance(icstate.particle.size());
	//fix chance
	MCVolumeAniso.setMinChance(0.1);
	MCVolumeAniso.setMaxChance(0.1);
#endif

	mover.addMove(MCSingleDispRot);
	mover.addMove(MCSingleRot);
	mover.addMove(MCSingleDisp);
	mover.addMove(MCSingleLocalDisp);
	mover.addMove(MCSingleLocalRot);
	mover.addMove(MCSingleNormal);
	mover.addMove(MCClusterMove);
#ifdef NPT
	mover.addMove(MCVolumeIso);
	mover.addMove(MCVolumeAniso);
#endif
	//s.checkSpace();
	int iterations = 10'000;
	for (int i = 0; i < iterations; i++)
	{
		double beta = betaMin + i*(betaMax-betaMin)/iterations;
		double pressure = pressureMin + i*(pressureMax-pressureMin) / iterations;

		s.setBeta(beta);
#ifdef NPT
		s.setPressure(pressure);
#endif
		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

		mover.execute(100'000);
		std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
		//s.checkSpace();
		cout << i << " "  << beta << " " << pressure << " " << s.getE() << " " << s.getV() << " ";
		for (int i = 0; i < mover.size(); i++)
		{
			cout << mover.getSuccessRate(i) << " ";
		}
		std::cout << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " ";
		std::cout << mover.getWork() << " ";
		cout << endl;
		vtf.save();
		mover.optimizeChances(1'000);
	}
}
