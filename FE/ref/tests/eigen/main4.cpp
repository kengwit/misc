#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
typedef Eigen::MatrixXd EgMatD2;
typedef Eigen::VectorXd EgVecD;

#include <blitz/array.h>
typedef blitz::Array<double,2> MatD2;
typedef blitz::Array<double,1> VecD;

	
	
int main()
{
	EgMatD2 Amat; Amat.resize(4,4);
	Amat << 1,2,3,4,
			5,6,7,8,
			9,10,11,12,
			13,14,15,16;
			
	cout << Amat.topLeftCorner(3,3) << endl;
	cout << "Amat.block(3,0,1,4) = " << endl;
	cout << Amat.block(3,0,1,4) << endl;
	Amat.block(3,0,1,4) << -1,-1,-1,-1;
	cout << "Amat.block(3,0,1,4) = " << endl;
	cout << Amat.block(3,0,1,4) << endl;
	
	EgVecD v1; v1.resize(4); 
	v1(0) = 2;
	v1(1) = 1;
	v1(2) = 1;
	v1(3) = 1.5;
	cout << v1.head(3).norm() << endl;
	cout << sqrt(6) << endl;
}
