#include <iostream>

#include "kdtree.h"


using namespace std;

bool print(int i)
{
	cout << i << endl;
}

int main()
{
	KdTree t;

	t.insert(0, {1,1,1});
	t.insert(1, {2,2,2});
	t.insert(2, {-1, -1, -1});
	t.insert(3, {-2,-2, -2});

	/*t.traverse(
			[](int n) {cout << n << endl; return (n != 2);}
			);*/

	/*t.traverse(
			[](int n) {cout << n << endl; return true;}, {-1,-3,-3}, {1,1,1}
			);*/

	t.remove(0, {1,1,1});
	/*t.traverse(
			[](int n) {cout << n << endl; return (true);}
			);*/
}
