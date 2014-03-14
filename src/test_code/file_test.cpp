#include <string>
#include <iostream>
#include <stdlib.h> // exit getenv

using namespace std;
 

int main() 
{
	char* data = getenv("DATA");
	string var = "P39900";
	string fname = data + string("validation_overlays/") + var + string(".sdf");
	cout << fname << endl;
	return 0;
}
