#include <iostream>

using namespace std;

int main() 
{
	int val = 5;
	int *point = &val;
	cout << val << endl;
	cout << "The pointer is: " << point << endl;
	cout << "The address of the val is: " << &val << endl;
	cout << "The dereferenced pointer is: " << *point << endl;
	return 0;
}
