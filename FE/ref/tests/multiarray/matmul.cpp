#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cassert>
#include <limits>
#include <random>
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
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0,1.0);
  
	Tensor2 A(boost::extents[3][3]);  
	Tensor2 B(boost::extents[3][3]);  
	Tensor2 C(boost::extents[3][3]);  
	Tensor2 CTranspose(boost::extents[3][3]);  
	Tensor2 BCt(boost::extents[3][3]);  
	Tensor2 ABC1(boost::extents[3][3]);  
	Tensor2 ABC2(boost::extents[3][3]);  
	
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			A[i][j] = distribution(generator);
			B[i][j] = distribution(generator);
			C[i][j] = distribution(generator);
	
		}
	}
	
	matlib_transpose(C, CTranspose);
	matlib_matmult(B, CTranspose, BCt);
	zeroT2(ABC1);
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			for ( int k = 0; k < 3; ++k ) {
				ABC1[i][j] += A[i][k]*BCt[k][j];
			}
		}
	}
	cout << "ABC1 = \n";
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			cout << ABC1[i][j] << " ";
		}
		cout << endl;
		
	}
	
	zeroT2(ABC2);
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			for ( int k = 0;  k < 3; ++k ) {
				for ( int l = 0; l < 3; ++l ) {
					ABC2[i][j] += A[i][k]*B[k][l]*C[j][l];
                }
			}
		}
	}
	cout << "ABC2 = \n";
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			cout << ABC2[i][j] << " ";
		}
		cout << endl;
		
	}
	
		
}

