#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cassert>
#include <limits>
using namespace std;
#include "boost/multi_array.hpp"

typedef boost::multi_array<double, 4> Tensor4;
typedef boost::multi_array<double, 2> Tensor2;
typedef boost::multi_array<double, 3> MatrixD3;
typedef boost::multi_array<double, 2> MatrixD2;
typedef boost::multi_array<double, 1> VectorD;

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


int main()
{
	int print_width = 30;
	int precision = 24;
	Tensor2 A1(boost::extents[3][3]);  
	Tensor2 A1T(boost::extents[3][3]);  
	
	double A2[3][3];
	double A2T[3][3];
	
	for ( int count = 0; count < 100000; ++count )
	{
	for ( int i = 0; i < 3; i++ ) {
		for ( int j = 0; j < 3; j++ ) {
			A1[i][j] = (double)(i)/(double)(j+sqrt(1.0/3.0));
		}
	}
	matlib_transpose(A1,A1T);
	
	for ( int i = 0; i < 3; i++ ) {
		for ( int j = 0; j < 3; j++ ) {
			A2[i][j] = (double)(i)/(double)(j+sqrt(1.0/3.0));
		}
	}
	for ( int i = 0; i < 3; i++ ) {
		for ( int j = 0; j < 3; j++ ) {
			A2T[i][j] = A2[j][i];
		}
	}
	
	cout << "A1 = " << endl;
	for ( int i = 0; i < 3; i++ ) {
		for ( int j = 0; j < 3; j++ ) {
			cout.width(print_width); 
			cout << std::fixed;
			cout.precision(precision);
			cout << std::right << A1[i][j] << " ";
		}
		cout << endl;
	}
	
	cout << "A1T = " << endl;
	for ( int i = 0; i < 3; i++ ) {
		for ( int j = 0; j < 3; j++ ) {
			cout.width(print_width); 
			cout << std::fixed;
			cout.precision(precision);
			cout << std::right << A1T[i][j] << " ";
		}
		cout << endl;
	}
	
	cout << "A2 = " << endl;
	for ( int i = 0; i < 3; i++ ) {
		for ( int j = 0; j < 3; j++ ) {
			cout.width(print_width); 
			cout << std::fixed;
			cout.precision(precision);
			cout << std::right << A2[i][j] << " ";
		}
		cout << endl;
	}
	cout << "A2T = " << endl;
	for ( int i = 0; i < 3; i++ ) {
		for ( int j = 0; j < 3; j++ ) {
			cout.width(print_width); 
			cout << std::fixed;
			cout.precision(precision);
			cout << std::right << A2T[i][j] << " ";
		}
		cout << endl;
	}
	
	double C2[3][3];
	for ( int i = 0; i < 3; i++ ) {
		for ( int j = 0; j < 3; j++ ) {
			C2[i][j] = 0.0;
		}
	}
	
	for ( int i = 0; i < 3; i++ ) {
		for ( int j = 0; j < 3; j++ ) {
			for ( int k = 0; k < 3; k++ ) {
				C2[i][j] = A1T[i][k]*A1[k][j];
			}
		}
	}
	
	cout << "C2 using A1\n";
	for ( int i = 0; i < 3; i++ ) {
		for ( int j = 0; j < 3; j++ ) {
			cout.width(print_width); 
			cout << std::fixed;
			cout.precision(precision);
			cout << std::right << C2[i][j] << " ";
		}
		cout << endl;
	}
	for ( int i = 0; i < 3; i++ ) {
		for ( int j = 0; j < 3; j++ ) {
			C2[i][j] = 0.0;
		}
	}
	for ( int i = 0; i < 3; i++ ) {
		for ( int j = 0; j < 3; j++ ) {
			for ( int k = 0; k < 3; k++ ) {
			C2[i][j] = A2T[i][k]*A2[k][j];
			}
		}
	}
	
	cout << "C2 using A2\n";
	for ( int i = 0; i < 3; i++ ) {
		for ( int j = 0; j < 3; j++ ) {
			cout.width(print_width); 
			cout << std::fixed;
			cout.precision(precision);
			cout << std::right << C2[i][j] << " ";
		}
		cout << endl;
	
	}
	}
	
}

