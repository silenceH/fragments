#include <iostream>

int main()
{
	int * a = new int;
	*a = 4;

	std::cout<<a<<"\t"<<*a<<std::endl;
	std::cout<<a<<"\t"<<&a<<std::endl;
}
