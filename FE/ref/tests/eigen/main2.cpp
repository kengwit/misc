#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <blitz/array.h>
typedef blitz::Array<double,2> MatD2;
typedef blitz::Array<double,1> VecD;

struct eigen_struct { 
	int index;
	double value;
	double vector[3];
	bool operator<(eigen_struct const &other) const { 
		return value < other.value;
	}
};
void solve_eigenvalue_3x3(MatD2 AMat, VecD BePr, MatD2 Nvec);

	
	
int main()
{
	MatD2 A(3,3);
	VecD  APr(3);
	MatD2 Nvec(3,3);
	
	A = 0;
	
	A = 1.2,1,0,
	     0,1,0,
		 0,0,1;
	
	cout << "A = " << endl;
	cout << A << endl;
	
	solve_eigenvalue_3x3(A, APr, Nvec);
	
	
	cout << "APr = " << endl;
	for ( int i = 0; i < 3; ++i ) {
		cout.precision(17);
		cout << APr(i) << ", vec = " << Nvec(i,0) << " " << Nvec(i,1) << " " << Nvec(i,2) << endl;
	}
	return 0;
}

void solve_eigenvalue_3x3(MatD2 AMat, VecD BePr, MatD2 Nvec)
{
	
	

	Eigen::EigenSolver<Eigen::Matrix3d> es;
	Eigen::Matrix3d A;
	
	// assign AMat to A
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			A(i,j) = AMat(i,j);
		}
	}
	
	es.compute(A);
	
	vector<eigen_struct> EigArray(3);
	
	//cout << "es.eigenvectors().col(0)[2] = " << endl << es.eigenvectors().col(0)[2] << endl;
	
	for ( int i = 0; i < 3; ++i )
	{
		EigArray[i].index = i;
		//cout << "es.eigenvalues()[i] = " << es.eigenvalues()[i].real() << endl;
		EigArray[i].value =  es.eigenvalues()[i].real();
		for ( int j = 0; j < 3; ++j )
		{
			EigArray[i].vector[j] = es.eigenvectors().col(i)[j].real();
		}
	}
	
	// sort based on ascending order of eigenvalues
	sort(EigArray.begin(),EigArray.end());
	
	// print using sorted indices
	//cout << "sorted eigenvalues: " << endl;
	//for ( int i = 0; i < 3; ++i )
	//{
	//	cout << EigArray[i].index << ", " << EigArray[i].value << ", eigevec: (" 
	//	     << EigArray[i].vector[0] << "," <<  EigArray[i].vector[1] << "," << EigArray[i].vector[2] << ")\n";
	//}
	
	// assign to output
	for ( int i = 0; i < 3; ++i )
	{
		BePr(i) = EigArray[i].value;
		for ( int j = 0; j < 3; ++j ) {
			Nvec(i,j) = EigArray[i].vector[j];
		}
	}
	
}

/*void solve_eigenvalue_3x3_selfadjoint(MatD2 AMat, VecD BePr, MatD2 Nvec)
{
	
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
	Eigen::Matrix3d A;
	
	// assign AMat to A
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			A(i,j) = AMat(i,j);
		}
	}
	
	es.compute(A);

	vector<eigen_struct> EigArray(3);
	
	//cout << "es.eigenvectors().col(0)[2] = " << endl << es.eigenvectors().col(0)[2] << endl;
	
	for ( int i = 0; i < 3; ++i )
	{
		EigArray[i].index = i;
		EigArray[i].value = es.eigenvalues()[i];
		for ( int j = 0; j < 3; ++j )
		{
			EigArray[i].vector[j] = es.eigenvectors().col(i)[j];
		}
	}
	
	// sort based on ascending order of eigenvalues
	sort(EigArray.begin(),EigArray.end());
	
	// print using sorted indices
	//cout << "sorted eigenvalues: " << endl;
	//for ( int i = 0; i < 3; ++i )
	//{
	//	cout << EigArray[i].index << ", " << EigArray[i].value << ", eigevec: (" 
	//	     << EigArray[i].vector[0] << "," <<  EigArray[i].vector[1] << "," << EigArray[i].vector[2] << ")\n";
	//}
	
	// assign to output
	for ( int i = 0; i < 3; ++i )
	{
		BePr(i) = EigArray[i].value;
		for ( int j = 0; j < 3; ++j ) {
			Nvec(i,j) = EigArray[i].vector[j];
		}
	}
	
}*/

