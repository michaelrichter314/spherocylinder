#include <string>
#include "analyzer.h"

using namespace fp;
using namespace std;

int main(int args, char** arg)
{
	Analyzer a( arg[1] );
	double theta = atof(arg[2]);
	a.setTheta(theta);

	do
	{
		cout << a.getAssociationSpherocyl() << " " << a.getAssociationEnd() << " " <<  a.getOrderParameter() << " " << a.getAvNN() << " " << a.getSmect() << " " << a.getSmecticParameter() << endl;
	} while (a.goToNextTimestep());

	return 0;
}
