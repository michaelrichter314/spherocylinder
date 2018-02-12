#include <iostream>
#include "cosinator.h"

using namespace std;

int main()
{
	Cosinator c;
	for (double x = -2; x <= 2; x+= 0.01)
	{
		cout << x << " " << c.acos(x) << endl;
	}
}
