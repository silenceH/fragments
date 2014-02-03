#include <iostream>
#include <math.h>

void getSinCos(double dx, double &dsin, double &dcos)
{
	dsin = sin(dx);
	dcos = cos(dx);
}

int main()
{
	double dsin = 0.0;
	double dcos = 0.0;
	
	getSinCos(3.14159, dsin, dcos);
	std::cout << "The sin value is: " << dsin << std::endl;
	std::cout << "The cos value is: " << dcos << std::endl;

	return 0;
}

