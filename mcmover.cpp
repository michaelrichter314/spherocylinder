#include "mcmover.h"

void McMove::setFunction(McMoveFunction f)
{
	function = f;
}

void McMove::setParameter(double p)
{
	parameters = 1;
	parameter1 = p;
}

void McMove::setParameter(double p1, double p2)
{
	parameters = 2;
	parameter2.p1 = p1;
	parameter2.p2 = p2;
}

void McMove::setSimulation(Simulation* s)
{
	sim = s;
}

void McMove::execute()
{
	switch (parameters)
	{
		case 0:
			successes += (*sim.*function)(NULL);
			break;
		case 1:
			successes += (*sim.*function)(&parameter1);
			break;
		case 2:
			successes += (*sim.*function)(&parameter2);
			break;
	}

	trials++;
}

void McMove::enable()
{
	enabled = true;
}

void McMove::disable()
{
	enabled = false;
}

bool McMove::isEnabled() const
{
	return enabled;
}

void McMove::setChance(double c)
{
	chance = c;
}

double McMove::getChance() const
{
	return chance;
}

McMove::McMove()
{

}

double McMove::getSuccessRate() const
{
	return (double)successes/trials;
}

void McMove::resetSuccessRate()
{
	trials = 0;
	successes = 0;
}

void McMove::setImportance(double i)
{
	importance = i;
}

double McMove::getImportance() const
{
	return importance;
}

void McMove::setMinChance(double m)
{
	minChance = m;
}

void McMove::setMaxChance(double m)
{
	maxChance = m;
}

double McMove::getMinChance() const
{
	return minChance;
}

double McMove::getMaxChance() const
{
	return maxChance;
}


// executes (iterations) number of monte carlo moves
// returns nothing
void McMover::execute(unsigned long iterations)
{
	if (needsToPrepare) {prepare();}

	for (unsigned long it = 0; it < iterations; it++)
	{
		double r = random.random(); //0..1

		// select a move
		int whichMove = activeMoves - 1;
		for (int i = 0; i < activeMoves-1; i++)
		{
			if (r <= maxRndValue[i])
			{
				whichMove = i;
				break;
			}
		}

		// execute the move
		activeMove[whichMove] -> execute();
	}
}

// adds new move to list
// returns id of move
int McMover::addMove(McMove m)
{
	move.push_back(m);
	needsToPrepare = true;
	return move.size() - 1;
}

// enables move with id (id)
// returns nothing
void McMover::enable(int id)
{
#ifdef debug
	if (id < 0 || id >= move.size())
	{
		std::cerr << "McMover::enable(" << id << "): id outside range." << std::endl;
		throw 1;
	}
#endif

	move[id].enable();
	needsToPrepare = true;
}

// disables move with id (id)
// returns nothing
void McMover::disable(int id)
{
#ifdef debug
	if (id < 0 || id >= move.size())
	{
		std::cerr << "McMover::disable(" << id << "): id outside range." << std::endl;
		throw 1;
	}
#endif

	move[id].disable();
	needsToPrepare = true;
}

// sets relative chance for move with (id)
// chance must be non-negative, but can be >=1
// returns nothing
void McMover::setChance(int id, double chance)
{
#ifdef debug
	if (id < 0 || id >= move.size())
	{
		std::cerr << "McMover::setChance(" << id << ", " << chance << "): id outside range." << std::endl;
		throw 1;
	}

	if (chance < 0.0)
	{
		std::cerr << "McMover::setChance(" << id << ", " << chance << "): chance outside range (0..+inf)." << std::endl;
		throw 1;
	}
#endif

	move[id].setChance(chance);
	needsToPrepare = true;
}

// gets percentage of successful moves of move (id)
// returns success rate
double McMover::getSuccessRate(int id)
{
#ifdef debug
	if (id < 0 || id >= move.size())
	{
		std::cerr << "McMover::getSuccessRate(" << id << "): id outside range." << std::endl;
		throw 1;
	}
#endif

	return move[id].getSuccessRate();
}

// resets success counter for move (id)
// returns nothing
void McMover::resetSuccessRate(int id)
{
#ifdef debug
	if (id < 0 || id >= move.size())
	{
		std::cerr << "McMover::resetSuccessRate(" << id << "): id outside range." << std::endl;
		throw 1;
	}
#endif

	move[id].resetSuccessRate();
}

// resets all success counters
// returns nothing
void McMover::resetSuccessRates()
{
	for (auto&& m : move)
	{
		m.resetSuccessRate();
	}
}

void McMover::prepare()
{
	if (!needsToPrepare) {return;}

	// look at active moves only
	activeMove.clear();

	for (auto&& m : move)
	{
		if (m.isEnabled())
		{
			activeMove.push_back(&m);
		}
	}
	activeMoves = activeMove.size();


	// sum up chances
	double sumOfChances = 0.0;
	for (const auto& m : activeMove)
	{
		sumOfChances += m->getChance();
	}

#ifdef debug
	if (activeMoves == 0 || sumOfChances <= 0.0)
	{
		//TODO: this is a problem!
	}
#endif

	// calculate maxRndValues
	maxRndValue.clear();
	maxRndValue.resize(activeMoves);
	maxRndValue[0] = activeMove[0]->getChance() / sumOfChances;

	for (int i = 1; i < activeMoves-1; i++)
	{
		maxRndValue[i] = maxRndValue[i-1] + activeMove[i]->getChance()/sumOfChances;
	}
	// fix rounding errors
	maxRndValue[activeMoves-1] = 1.0;

	// done! time for a beer
	needsToPrepare = false;
}

// side effect: resets success rates
void McMover::optimizeChances(unsigned long iterations)
{
	resetSuccessRates();

	std::vector<double> runtime(move.size());

	// test all mc step types
	for (int i = 0; i < move.size(); i++)
	{
		if (move[i].getMinChance() < move[i].getMaxChance() && move[i].isEnabled())
		{
			auto start = std::chrono::high_resolution_clock::now();
			for (unsigned long it = 0; it < iterations; it++)
			{
				move[i].execute();
			}
			auto end = std::chrono::high_resolution_clock::now();

			runtime[i] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
		}
		else
		{
			runtime[i] = 1.0;
		}
	}

	std::vector<double> workDone(move.size());
	std::vector<double> newWeight(move.size());
	for (int i = 0; i < workDone.size(); i++)
	{
		workDone[i] = move[i].getImportance() * move[i].getSuccessRate();
		newWeight[i] = workDone[i] / runtime[i];
	}
	double maxWeight = 0;
	for (int i = 0; i < newWeight.size(); i++)
	{
		if (newWeight[i] > maxWeight && move[i].getMinChance() < move[i].getMaxChance())
		{
			maxWeight = newWeight[i];
		}
	}
	for (int i = 0; i < newWeight.size(); i++)
	{
		newWeight[i] *= 100.0 / maxWeight;
		newWeight[i] = std::fmax(newWeight[i], move[i].getMinChance());
		newWeight[i] = std::fmin(newWeight[i], move[i].getMaxChance());
		move[i].setChance( ( newWeight[i] > move[i].getChance() ? newWeight[i] : 0.25*newWeight[i]+0.75*move[i].getChance() ) ); // increasing is fast, decreasing slow
	}

	resetSuccessRates();
}

fp::Randomizer* McMover::getRandomizer()
{
	return &random;
}

int McMover::size() const
{
	return move.size();
}

double McMover::getWork() const
{
	double toReturn = 0.0;
	for (int i = 0; i < move.size(); i++)
	{
		if (move[i].isEnabled())
		{
			toReturn += move[i].getSuccessRate()*move[i].getImportance()*move[i].getChance();
		}
	}

	return toReturn;
}

double McMover::getChance(int i) const
{
	return move[i].getChance();
}







