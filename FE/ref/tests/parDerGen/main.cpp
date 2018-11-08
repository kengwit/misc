#include <iostream>
#include <iomanip>
using namespace std;

#include "ParDerGen.h"

int main()
{
	Eigen::Matrix3d X; 
	//X << 1.3,1,1,
	//    1,4.5,1,
	//	 1,1,3.3;
	
	X << 1,0,0,
		 0,0.99,0,
		 0,0,1.03;
		 
	cout << "X:" << endl << X << endl;
	
	Eigen::MatrixXd L;
	ParDerGen(X, L);
	cout << std::setprecision(5);
	cout << std::scientific;
	cout << "L: " << endl  << L << endl;
	
	

}