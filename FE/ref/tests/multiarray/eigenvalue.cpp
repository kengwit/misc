#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cassert>
#include <limits>
#include <random>
#include <chrono>
using namespace std;

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
typedef Eigen::MatrixXd EgMatD2;
typedef Eigen::VectorXd EgVecD;

#include "boost/multi_array.hpp"
typedef boost::multi_array<double, 4> Tensor4;
typedef boost::multi_array<double, 2> Tensor2;
typedef boost::multi_array<double, 3> MatrixD3;
typedef boost::multi_array<double, 2> MatrixD2;
typedef boost::multi_array<double, 1> VectorD;

struct eigen_struct { 
	int index;
	double value;
	double vector[3];
	bool operator<(eigen_struct const &other) const { 
		return value < other.value;
	}
};

void zeroT2(Tensor2 & A)
{
	fill(A.origin(), A.origin()+A.num_elements(), 0.0);
}

void matlib_transpose(Tensor2 & A, Tensor2 & AT)
{
	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
			AT[i][j] = A[j][i];
}

void matlib_matmult(Tensor2 & A, Tensor2 & B, Tensor2 & C)
{
	zeroT2(C);
	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
			for ( int k = 0; k < 3; ++k )
				C[i][j] += A[i][k]*B[k][j];
}


void solve_eigenvalue_3x3(Tensor2 & AMat, VectorD & APr, MatrixD2 & Nvec)
{
	
	Eigen::EigenSolver<Eigen::Matrix3d> es;
	Eigen::Matrix3d A;
	
	// assign AMat to A
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			A(i,j) = AMat[i][j];
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
		
		//cout << "eigenvec " << i << " = " << endl;
		//cout << es.eigenvectors().col(i) << endl;
		//cout << es.eigenvectors().col(i).norm() << endl;
		
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
		APr[i] = EigArray[i].value;
		for ( int j = 0; j < 3; ++j ) {
			Nvec[i][j] = EigArray[i].vector[j];
		}
		
		/*double norm = 0.0;
		for ( int j = 0; j < 3; ++j ) {
			cout << Nvec[i][j] << " ";
			norm += Nvec[i][j]*Nvec[i][j];
		}
		norm = sqrt(norm);
		cout << "norm " << i << " = " << norm << endl;*/
	}
	
	/**/
	
}


int main()
{
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator (seed);
	std::uniform_real_distribution<double> distribution(0.0,1.0);
  
	Tensor2 A(boost::extents[3][3]);  
	Tensor2 Asym(boost::extents[3][3]);  
	VectorD  AsymPr(boost::extents[3]);      // eigenvalues of Be_TR
	MatrixD2 Nvec_Asym(boost::extents[3][3]);   // eigenvectors of Be_TR
	MatrixD3 M_Asym(boost::extents[3][3][3]);   // eigenvectors of Be_TR
	Tensor2 AsymRecon(boost::extents[3][3]);  
	
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			A[i][j] = (double)(i+j);//distribution(generator);
		}
	}
	
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			Asym[i][j] = 0.5*(A[i][j]+A[j][i]);
		}
	}
	
	cout << "Asym = \n";
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			cout << Asym[i][j] << " ";
		}
		cout << endl;
	}
	solve_eigenvalue_3x3(Asym, AsymPr, Nvec_Asym);
	
	cout << "eigenvectors are = " << endl;
	for ( int i = 0; i < 3; ++i ) {
		cout << "Nvec " << i << " = ";
		for ( int j = 0; j < 3; ++j ) {
			cout << Nvec_Asym[i][j] << " ";
		}		
		cout << endl;
	}

	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			M_Asym[0][i][j] = Nvec_Asym[0][i]*Nvec_Asym[0][j];
			M_Asym[1][i][j] = Nvec_Asym[1][i]*Nvec_Asym[1][j];
			M_Asym[2][i][j] = Nvec_Asym[2][i]*Nvec_Asym[2][j];
		}		
	}
	
	
	zeroT2(AsymRecon);
	for ( int k = 0; k < 3; ++k ) {
		for ( int i = 0; i < 3; ++i ) {
			for ( int j = 0; j < 3; ++j ) {
				AsymRecon[i][j] += AsymPr[k]*Nvec_Asym[k][i]*Nvec_Asym[k][j];
			}		
		}
	}
	
	/*zeroT2(AsymRecon);
	for ( int k = 0; k < 3; ++k ) {
		for ( int i = 0; i < 3; ++i ) {
			for ( int j = 0; j < 3; ++j ) {
				AsymRecon[i][j] += AsymPr[k]*M_Asym[k][i][j];
			}		
		}
	}*/
	
	cout << "AsymRecon = \n";
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			cout << AsymRecon[i][j] << " ";
		}
		cout << endl;
	}
	
	cout << "zero AsymRecon\n";
	//AsymRecon = 0.0;
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			cout << AsymRecon[i][j] << " ";
		}
		cout << endl;
		
	}
	
}

