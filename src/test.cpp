
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h> // exit getenv

using namespace std;

int main()
{
	char* data = getenv("DATA");
	string fname = data + string("validation_overlays/") +"P39900" + string(".sdf");
	cout << fname << endl;


}
