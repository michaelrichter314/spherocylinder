#pragma once
#include "fpspace.h"
#include <fpdim3.h>
#include <cmath>
#include "fpbinarysearchtree.h"
#ifdef debug
#include <iostream>
#endif

namespace fp
{

	// space size and grid size
	struct GridSize
	{
		int3 gridSize;
		double3 spaceSize;
	};

	class Grid : public Space
	{
		public:
			Grid();
			Grid(int3, double, double3);
			void initialize(int3, double, double3);

			// stuff to extract occupancy information in one grid cell + neighboring grid cells
			void getIDsNear(const double3&);
			int getNextID();
			bool gotMoreIDs();

			// for saving a copy of the grid
			int3 getMaxCells() const;

			// ovverride stuff
			void changeVolume(double3);
			void insert(Spherocyl&);
			void remove(Spherocyl&);
//			void move(int, const double3&, const double3&);
		private:
			int3 maxCells;
			double minMetric;

			// for occupancy info stuff
			BinarySearchTree<int>* currentCell = NULL;
			int3 currentCellOrigin;
			int3 currentCellD;
			int3 gridCellPBC(int3) const;
			void goToNextCell();
			bool currentCellIsLastCell = false;
			int3 minCellD;
			int3 maxCellD;
			

			vector3<BinarySearchTree<int>> tree;

			double3 cellSize;
			int3 cells;

			bool isFullX;
			bool isFullY;
			bool isFullZ;

			bool isFull(double, int) const;

			int3 getGridCell(double3) const;
	};
}

