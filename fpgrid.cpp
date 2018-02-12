#include "fpgrid.h"

namespace fp
{
	Grid::Grid()
	{
	}


	Grid::Grid(int3 maxcells, double minmetric, double3 metricsize)
	{
		initialize(maxcells, minmetric, metricsize);
	}


	void Grid::initialize(int3 maxcells, double minmetric, double3 metricsize)
	{
		maxCells = maxcells;
		minMetric = minmetric;
		spaceSize = metricsize;

		tree.clear();
		tree.resize(maxCells);
		changeVolume(spaceSize);
	}


	void Grid::changeVolume(double3 newVolume)
	{
		// determine type of tree for each dimension
		isFullX = isFull(newVolume.x, maxCells.x);
		isFullY = isFull(newVolume.y, maxCells.y);
		isFullZ = isFull(newVolume.z, maxCells.z);

		// calculate cells in each direction
		cells.x = ( isFullX ? maxCells.x : (int)(newVolume.x/minMetric) );
		cells.y = ( isFullY ? maxCells.y : (int)(newVolume.y/minMetric) );
		cells.z = ( isFullZ ? maxCells.z : (int)(newVolume.z/minMetric) );

		// calculate tree size for each dimension
		cellSize.x = ( isFullX ? newVolume.x/cells.x : minMetric );
		cellSize.y = ( isFullY ? newVolume.y/cells.y : minMetric );
		cellSize.z = ( isFullZ ? newVolume.z/cells.z : minMetric );

		spaceSize = newVolume;
	}


	void Grid::insert(Spherocyl& s)
	{
		int ID = s.getID();
		const double3& p = s.getP();
		tree.at( getGridCell(p) ).insert(ID);
	}


	void Grid::remove(Spherocyl& s)
	{
		int ID = s.getID();
		const double3& p = s.getP();
		auto& v = tree.at( getGridCell(p) );

		v.remove(ID);
	}

/*	void Grid::move(int ID, const double3& oldP, const double3& newP)
	{
		int3 oldG = getGridCell(oldP);
		int3 newG = getGridCell(newP);

		// only move if tree cell changed
		if (oldG != newG)
		{
			tree.at(oldG).remove(ID);
			tree.at(newG).insert(ID);
		}
	}*/


	bool Grid::isFull(double size, int cells) const
	{
		return ( size > cells*minMetric );
	}


	int3 Grid::getGridCell(double3 p) const
	{
		int3 toReturn;

		toReturn.x = ( isFullX ? (int)(cells.x * p.x / spaceSize.x) : std::min((int)(p.x / minMetric), cells.x-1) );
		toReturn.y = ( isFullY ? (int)(cells.y * p.y / spaceSize.y) : std::min((int)(p.y / minMetric), cells.y-1) );
		toReturn.z = ( isFullZ ? (int)(cells.z * p.z / spaceSize.z) : std::min((int)(p.z / minMetric), cells.z-1) );

#ifdef debug
		if (p.x >= spaceSize.x || p.x < 0.0) { std::cerr << "getGridCell(): invalid position " << p << std::endl; }
		if (p.y >= spaceSize.y || p.y < 0.0) { std::cerr << "getGridCell(): invalid position " << p << std::endl; }
		if (p.z >= spaceSize.z || p.z < 0.0) { std::cerr << "getGridCell(): invalid position " << p << std::endl; }

		if (toReturn.x >= cells.x) { std::cerr << "getGridCell(): invalid tree " << toReturn << ", should be < " << cells << std::endl; }
		if (toReturn.y >= cells.y) { std::cerr << "getGridCell(): invalid tree " << toReturn << ", should be < " << cells << std::endl; }
		if (toReturn.z >= cells.z) { std::cerr << "getGridCell(): invalid tree " << toReturn << ", should be < " << cells << " with minMetric=" << minMetric << " p.z = " << p.z << std::endl; }

#endif

		return toReturn;
	}




	void Grid::getIDsNear(const double3& p)
	{
		currentCellOrigin = getGridCell(p);
		minCellD = {(cells.x > 1 ? -1 : 0), (cells.y > 1 ? -1 : 0), (cells.z > 1 ? -1 : 0)};
		maxCellD = {(cells.x > 2 ? 1 : 0), (cells.y > 2 ? 1 : 0), (cells.z > 2 ? 1 : 0)};
		currentCellD = minCellD;

		currentCellIsLastCell = (cells.x == 1 && cells.y == 1 && cells.z == 1);
		
		//currentCell = &(tree.at( {0,0,0} ));
		currentCell = &(tree.at(gridCellPBC(currentCellOrigin + currentCellD) ));
		currentCell -> startTraversal();
		while (!(currentCell -> gotAnotherElement()) && !currentCellIsLastCell)
		{
			goToNextCell();
			currentCell = &(tree.at(gridCellPBC(currentCellOrigin + currentCellD) ));
			currentCell -> startTraversal();
		}
	}

	int Grid::getNextID()
	{
		int toReturn = currentCell -> getNextElement();

#ifdef debug
		std::cout << currentCellOrigin << std::endl;
		std::cout << currentCellD<< std::endl;
		std::cout << cells << std::endl;
		std::cout << currentCellOrigin+currentCellD << std::endl;
		std::cout << gridCellPBC(currentCellOrigin+currentCellD) << std::endl;
#endif
		while (!(currentCell -> gotAnotherElement()) && !currentCellIsLastCell)
		{
			goToNextCell();
			currentCell = &(tree.at(gridCellPBC(currentCellOrigin + currentCellD) ));
			currentCell -> startTraversal();
		}



		return toReturn;
	}

	void Grid::goToNextCell()
	{
		if (currentCellD.x == maxCellD.x)
		{
			currentCellD.x = minCellD.x;
			currentCellD.y++;

			if (currentCellD.y == maxCellD.y + 1)
			{
				currentCellD.y = minCellD.y;
				currentCellD.z++;

			}
		}
		else
		{
			currentCellD.x++;
			if (currentCellD == maxCellD) {currentCellIsLastCell = true;}
		}
	}


	bool Grid::gotMoreIDs()
	{
		return currentCell -> gotAnotherElement();
	}	

	int3 Grid::getMaxCells() const
	{
		return maxCells;
	}

	int3 Grid::gridCellPBC(int3 in) const
	{
		if (in.x < 0) {in.x += cells.x * (-in.x/cells.x + 1);}
		if (in.y < 0) {in.y += cells.y * (-in.y/cells.y + 1);}
		if (in.z < 0) {in.z += cells.z * (-in.z/cells.z + 1);}

		in.x = in.x % cells.x;
		in.y = in.y % cells.y;
		in.z = in.z % cells.z;

		return in;

	}
}
