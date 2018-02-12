#include "kdtree.h"

void KdTree::initialize(const double3& size, Simulation* s)
{
	clear();
	spaceSize = size;
	sim = s;
}

void KdTree::insert(int N, const double3& p)
{
	
	// create point
	KdPoint* point = new KdPoint;
	point->N = N;
	point->p = p;

	// insert in correct position
	if (root == NULL)
	{
		root = point;
		point -> level = 0;
		return;
	}

	// traverse tree
	int level = 0;
	KdPoint* pos = root;
	bool done = false;
	bool goesOnLeft;

	while (!done)
	{
		switch (level%3)
		{
			case 0: //x-dimension
				goesOnLeft = (p.x < pos->p.x);
				break;
			case 1: //y-dimension
				goesOnLeft = (p.y < pos->p.y);
				break;
			case 2: //z-dimension
				goesOnLeft = (p.z < pos->p.z);
				break;
		}

		if (goesOnLeft)
		{
			// goes on left
			if (pos -> left != NULL)
			{
				pos = pos -> left;
			}
			else
			{
				// found it's place!
				pos->left = point;
				point -> parent = pos;
				point -> level = (level+1) % 3;
				done = true;
			}
		}
		else
		{
			// goes on right
			if (pos -> right != NULL)
			{
				pos = pos->right;
			}
			else
			{
				// found it's place!
				pos->right = point;
				point -> parent = pos;
				point -> level = (level+1) % 3;
				done = true;
			}

		}

		level++;
	}

}

void KdTree::remove(int N, const double3& p)
{
	remove(find(N,p));
}

void KdTree::remove(KdPoint* point)
{
	if (point == NULL)
	{
		return;
	}

	if (point -> left == NULL && point -> right == NULL)
	{
		// no children, so we can just delete this one
		// remove from parent, if there is one
		if (point -> parent == NULL)
		{
			delete point;
			root = NULL;
			return;
		}
		else
		{
			if (point -> parent -> left == point)
			{
				point -> parent -> left = NULL;
			}
			else
			{
				point -> parent -> right = NULL;
			}

			// delete it
			delete point;
		}
	}
	else if (point -> right != NULL)
	{
		auto replacement = findMin(point -> right, point -> level);
		point -> N = replacement -> N;
		point -> p = replacement -> p;

		remove(replacement); // is lower in the tree, so recursion will terminate
	}
	else
	{
		auto replacement = findMin(point -> left, point -> level);
		point -> N = replacement -> N;
		point -> p = replacement -> p;

		remove(replacement);

		// what was left is now right. works since there is no right subtree
		point -> right = point -> left;
		point -> left = NULL;

	}
}

KdPoint* KdTree::findMin(KdPoint* point, int level)
{
	if (point == NULL) {return NULL;}

	if (point -> level == level)
	{
		// only need to check left subtree
		auto left = findMin(point -> left, level);
		if (left != NULL)
		{
			return left;
		}
		else
		{
			return point;
		}
	}
	else
	{
		// need to check both sides. all is possible
		double myValue = getValue(point, level);
		auto left = findMin(point -> left, level);
		auto right = findMin(point -> right, level);

		auto currentBestPoint = point;
		auto currentBestValue = myValue;

		if (left != NULL && getValue(left, level) < currentBestValue)
		{
			currentBestPoint = left;
			currentBestValue = getValue(left, level);
		}

		if (right != NULL && getValue(right, level) < currentBestValue)
		{
			currentBestPoint = right;
		}

		return currentBestPoint;
	}
}

inline double KdTree::getValue(KdPoint* p, int l)
{
	return ((&(p->p.x))[l]);
	
	/*switch (l)
	{
		case 0:
			return p->p.x;
			break;
		case 1:
			return p->p.y;
			break;
		case 2:
			return p->p.z;
	}

	return 0.0;*/
}
inline double KdTree::getValue(const double3& p, int l)
{
	return ((&p.x)[l]);

	/*
	double r;
	switch (l)
	{
		case 0:
			r = p.x;
			break;
		case 1:
			r = p.y;
			break;
		case 2:
			r = p.z;
	}
	return r*/
}

KdPoint* KdTree::find(int N, const double3& p)
{
	if (root == NULL) {return NULL;}

	if (root->N == N) {return root;}

	// traverse tree. it's not root
	int level = 0;
	KdPoint* pos = root;
	bool goesOnLeft;

	while (true)
	{
		if (pos -> N == N) {return pos;}

		switch (level%3)
		{
			case 0: //x-dimension
				goesOnLeft = (p.x < pos->p.x);
				break;
			case 1: //y-dimension
				goesOnLeft = (p.y < pos->p.y);
				break;
			case 2: //z-dimension
				goesOnLeft = (p.z < pos->p.z);
				break;
		}

		if (goesOnLeft)
		{
			// goes on left
			if (pos -> left != NULL)
			{
				pos = pos->left;
			}
			else
			{
				// found it's place, but it's not there!
				std::cerr << "nothing's there! N=" << N << " " << p << std::endl;
				auto found = brute_force_find(N);
				if (found == NULL)
				{
					std::cerr << "not in tree!" << std::endl;
				}
				else
				{
					std::cerr << "found it here: " << found -> p << std::endl;
				}
				return found;
			}
		}
		else
		{
			// goes on right
			if (pos -> right != NULL)
			{
				pos = pos->right;
			}
			else
			{
				// found it's place, but it's not there!
				std::cerr << "nothing's there! N=" << N << " " << p << std::endl;
				auto found = brute_force_find(N);
				if (found == NULL)
				{
					std::cerr << "not in tree!" << std::endl;
				}
				else
				{
					std::cerr << "found it here: " << found -> p << std::endl;
				}
				return found;
			}

		}

		level++;
	}
}

KdPoint* KdTree::brute_force_find(int N)
{
	return brute_force_find(N, root);
}
	
KdPoint* KdTree::brute_force_find(int N, KdPoint* pos)
{
	if (pos == NULL) {return NULL;}
	if (pos -> N == N) {return pos;}
	auto l = brute_force_find(N, pos->left);
	if (l != NULL) {return l;}
	return brute_force_find(N, pos->right);
}
	

KdTree::~KdTree()
{
	clear();
}

void KdTree::clear()
{
	removeSubtree(root);
}

void KdTree::removeSubtree(KdPoint* p)
{
	if (p == NULL) {return;}

	removeSubtree(p->left);
	removeSubtree(p->right);

	if (root == p)
	{
		root = NULL;
	}

	delete(p);
}

void KdTree::traverse(TraverseFunction f)
{
	traverse(root, f);
}

bool KdTree::traverse(KdPoint* p, TraverseFunction f)
{
	if (p == NULL) {return true;}

	if (!traverse(p -> left, f)) {return false;}
	if (!(*sim.*f)(p->N)) {return false;}
	if (!traverse(p -> right, f)) {return false;}

	return true;
}


vector<int> KdTree::getInBox(const double3& ll, const double3& ur)
{
	vector<int> toReturn;

	// PBC stuff
	bool invX = ( ll.x > ur.x );
	bool invY = ( ll.y > ur.y );
	bool invZ = ( ll.z > ur.z );

	bool x,y,z;
	double3 l;
	double3 r;


	x = true;
	do
	{
		y = true;
		do
		{
			z = true;
			do
			{
				l.x = ( invX ? ( x ? 0.0 : ll.x ) : ll.x );
				l.y = ( invY ? ( y ? 0.0 : ll.y ) : ll.y );
				l.z = ( invZ ? ( z ? 0.0 : ll.z ) : ll.z );

				r.x = ( invX ? ( x ? ur.x : spaceSize.x ) : ur.x );
				r.y = ( invY ? ( y ? ur.y : spaceSize.y ) : ur.y );
				r.z = ( invZ ? ( z ? ur.z : spaceSize.z ) : ur.z );

				getInBox(root, l, r, toReturn);
				z=!z;
			} while ( z != invZ );

			y=!y;
		} while ( y != invY );
		x=!x;
	} while ( x != invX );

	return toReturn;
}

bool KdTree::traverse(TraverseFunction f, const double3& ll, const double3& ur)
{
	// PBC stuff
	bool invX = ( ll.x > ur.x );
	bool invY = ( ll.y > ur.y );
	bool invZ = ( ll.z > ur.z );

	bool x,y,z;
	double3 l;
	double3 r;

	x = true;
	do
	{
		y = true;
		do
		{
			z = true;
			do
			{
				l.x = ( invX ? ( x ? 0.0 : ll.x ) : ll.x );
				l.y = ( invY ? ( y ? 0.0 : ll.y ) : ll.y );
				l.z = ( invZ ? ( z ? 0.0 : ll.z ) : ll.z );

				r.x = ( invX ? ( x ? ur.x : spaceSize.x ) : ur.x );
				r.y = ( invY ? ( y ? ur.y : spaceSize.y ) : ur.y );
				r.z = ( invZ ? ( z ? ur.z : spaceSize.z ) : ur.z );

				if (!traverse(root, f, l, r)) {return false;}
				z=!z;
			} while ( z != invZ );

			y=!y;
		} while ( y != invY );
		x=!x;
	} while ( x != invX );

	return true;
}

bool KdTree::traverse(KdPoint* p, TraverseFunction f, const double3& ll, const double3& ur)
{
	if (p == NULL) {return true;}

	if (getValue(p, p->level) >= getValue(ll, p->level))
	{
		if (!traverse(p -> left, f, ll, ur)) {return false;}
	}

	if (p->p.x >= ll.x && p->p.y >= ll.y && p->p.z >= ll.z && p->p.x <= ur.x && p->p.y <= ur.y && p->p.z <= ur.z)
	{
		if (!(*sim.*f)(p->N)) {return false;}
	}

	if (getValue(p, p->level) <= getValue(ur, p->level))
	{
		if (!traverse(p -> right, f, ll, ur)) {return false;}
	}

	return true;
}

void KdTree::getInBox(KdPoint* p, const double3& ll, const double3& ur, vector<int>& toReturn)
{
	if (p == NULL) {return;}

	if (getValue(p, p->level) >= getValue(ll, p->level))
	{
		getInBox(p -> left, ll, ur, toReturn);
	}

	if (p->p.x >= ll.x && p->p.y >= ll.y && p->p.z >= ll.z && p->p.x <= ur.x && p->p.y <= ur.y && p->p.z <= ur.z)
	{
		toReturn.push_back(p->N);
	}

	if (getValue(p, p->level) <= getValue(ur, p->level))
	{
		getInBox(p -> right, ll, ur, toReturn);
	}
}

void KdTree::changeVolume(double3 newSize)
{
	spaceSize = newSize;
}

