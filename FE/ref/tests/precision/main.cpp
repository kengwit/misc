#include <iostream>
#include <cmath>
using namespace std;

int main()
{
	//double pt = sqrt(1.0/3.0);
    long double pt = sqrt(static_cast<long double>(1.0)/static_cast<long double>(3.0));
		
	cout << "pt = ";
	cout << std::fixed;
	cout.precision(30);
	cout << pt << endl;
	
	return 0;
}

