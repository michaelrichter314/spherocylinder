#pragma once

#include <iostream>
#include <vector>
#include <chrono>
#include <fprandom.h>

class Simulation;

using McMoveFunction = bool (Simulation::*)(const void*);

struct McMoveParameter2
{
	double p1;
	double p2;
};

class McMove
{
	protected:
		McMoveFunction function;
		bool enabled = false;
		double chance = 0.0;

		unsigned long trials = 0;
		unsigned long successes = 0;
		int parameters = 0;
		double parameter1;
		McMoveParameter2 parameter2;
		Simulation* sim;

		double importance = 1.0;
		double minChance = 1.0;
		double maxChance = 100.0;

	public:
		void setFunction(McMoveFunction);
		void setParameter(double);
		void setParameter(double, double);
		void setSimulation(Simulation*);

		void execute();

		void enable();
		void disable();
		bool isEnabled() const;

		void setChance(double);
		double getChance() const;

		McMove();

		double getSuccessRate() const;
		void resetSuccessRate();

		void setImportance(double);
		double getImportance() const;

		double getMinChance() const;
		double getMaxChance() const;
		void setMinChance(double);
		void setMaxChance(double);
};


class McMover
{
	public:
		int addMove(McMove);
		void enable(int);
		void disable(int);
		void setChance(int, double);

		void execute(unsigned long);

		double getSuccessRate(int);
		double getChance(int) const;
		void resetSuccessRate(int);
		void resetSuccessRates();

		void optimizeChances(unsigned long);

		//TODO: void optimizeParameters(unsigned long);

		fp::Randomizer* getRandomizer();

		int size() const;
		double getWork() const;

	protected:
		std::vector<McMove> move;

		// vector to decide which move is executed, contains maximal number random.random() can return which
		// leads to the move being executed
		// only contains active moves
		std::vector<double> maxRndValue;
		int activeMoves = 0;
		std::vector<McMove*> activeMove;
		
		// if the mover has been changed, it needs to prepare before executing
		// is called by execute()
		void prepare();
		bool needsToPrepare = true;

		fp::Randomizer random;
};
