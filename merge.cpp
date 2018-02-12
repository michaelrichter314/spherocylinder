#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int args, char** arg)
{
	if (args < 3)
	{
		cerr << "usage: " << arg[0] << " <file1> <file2>" << endl;
		return 1;
	}

	ifstream file1(arg[1]);
	ifstream file2(arg[2]);
	
	string buffer1;
	string buffer2;


	while (getline(file1, buffer1) && getline(file2, buffer2))
	{
		cout << buffer1 << " " << buffer2 << endl;
	}

	return 0;
}
