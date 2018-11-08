#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cassert>
#include <limits>
#include <random>
#include <chrono>
using namespace std;

#include "boost/multi_array.hpp"
typedef boost::multi_array<double, 4> Tensor4;
typedef boost::multi_array<double, 2> Tensor2;
typedef boost::multi_array<double, 3> MatrixD3;
typedef boost::multi_array<double, 2> MatrixD2;
typedef boost::multi_array<double, 1> VectorD;

template<typename T>
void zero2(T & A)
{
	fill(A.origin(), A.origin()+A.num_elements(), 0.0);
}

template<typename T>
void zero4(T & A);
{
	fill(A.origin(), A.origin()+A.num_elements(), 0.0);
}

void matlib_transpose(Tensor2 & A, Tensor2 & AT)
{
	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
			AT[i][j] = A[j][i];
}

void matlib_matmult9(MatrixD2 & A, MatrixD2 & B, MatrixD2 & C)
{
	zero2(C);
	for ( int i = 0; i < 9; ++i )
		for ( int j = 0; j < 9; ++j )
			for ( int k = 0; k < 9; ++k )
				C[i][j] += A[i][k]*B[k][j];
}

void converTensor4ToMat(Tensor4 & T, MatrixD2 & A);

int main()
{
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator (seed);
	std::uniform_real_distribution<double> distribution(0.0,1.0);
  
	Tensor4 A(boost::extents[3][3][3][3]);  
	Tensor4 B(boost::extents[3][3][3][3]);  
	Tensor4 C(boost::extents[3][3][3][3]);  
	
	MatrixD2 AMat(boost::extents[9][9]);  
	MatrixD2 BMat(boost::extents[9][9]);  
	MatrixD2 CMat(boost::extents[9][9]);  
	
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			for ( int k = 0; k < 3; ++k ) {
				for ( int l = 0; l < 3; ++l ) {
					A[i][j][k][l] = distribution(generator);
					B[i][j][k][l] = distribution(generator);
					C[i][j][k][l] = distribution(generator);
				}
			}
		}
	}
	
	converTensor4ToMat(A, AMat);
	converTensor4ToMat(B, BMat);
	converTensor4ToMat(C, CMat);
	
	// tensor product AB
	Tensor4 ABC(boost::extents[3][3][3][3]);  
	zero4(ABC);
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			for ( int k = 0; k < 3; ++k ) {
				for ( int l = 0; l < 3; ++l ) {
					for ( int p = 0; p < 3; ++p ) {
						for ( int q = 0; q < 3; ++q ) {
							for ( int m = 0; m < 3; ++m ) {
								for ( int n = 0; n < 3; ++n ) {
									ABC[i][j][k][l] += A[i][j][p][q]*B[p][q][m][n]*C[m][n][k][l];
								}
							}
						}
					}
				}
			}
		}
	}
	
	
	// matrix product AB
	MatrixD2 ABCMatFromTensor(boost::extents[9][9]);  
	converTensor4ToMat(ABC,ABCMatFromTensor);
	cout << "from tensor\n";
	for ( int i = 0; i < 9; ++i ) {
		for ( int j = 0; j < 9; ++j ) {
			cout << ABCMatFromTensor[i][j] << " ";
		}
		cout << endl;
	}
	
	// matrix product AB
	MatrixD2 ABMat(boost::extents[9][9]);  
	MatrixD2 ABCMat(boost::extents[9][9]);  
	matlib_matmult9(AMat, BMat, ABMat);
	matlib_matmult9(ABMat, CMat, ABCMat);
	cout << "from matrix\n";
	for ( int i = 0; i < 9; ++i ) {
		for ( int j = 0; j < 9; ++j ) {
			cout << ABCMat[i][j] << " ";
		}
		cout << endl;
	}
	
	MatrixD2 DiffMat(boost::extents[9][9]);  
	zero2(DiffMat);
	for ( int i = 0; i < 9; ++i ) {
		for ( int j = 0; j < 9; ++j ) {
			DiffMat[i][j] = ABCMat[i][j] - ABCMatFromTensor[i][j];
		}
	}
	
	double sum = 0.0;
	for ( int i = 0; i < 9; ++i ) {
		for ( int j = 0; j < 9; ++j ) {
			sum += DiffMat[i][j];
		}
	}
	
	cout << "sum = " << sum << endl;
	
}

void converTensor4ToMat(Tensor4 & T, MatrixD2 & A)
{
	A[0][0] = T[0][0][0][0];
	A[0][1] = T[0][0][1][1];
	A[0][2] = T[0][0][2][2];
	A[0][3] = T[0][0][0][1];
	A[0][4] = T[0][0][1][0];
	A[0][5] = T[0][0][1][2];
	A[0][6] = T[0][0][2][1];
	A[0][7] = T[0][0][2][0];
	A[0][8] = T[0][0][0][2];
	
	A[1][0] = T[1][1][0][0];
	A[1][1] = T[1][1][1][1];
	A[1][2] = T[1][1][2][2];
	A[1][3] = T[1][1][0][1];
	A[1][4] = T[1][1][1][0];
	A[1][5] = T[1][1][1][2];
	A[1][6] = T[1][1][2][1];
	A[1][7] = T[1][1][2][0];
	A[1][8] = T[1][1][0][2];
	
	A[2][0] = T[2][2][0][0];
	A[2][1] = T[2][2][1][1];
	A[2][2] = T[2][2][2][2];
	A[2][3] = T[2][2][0][1];
	A[2][4] = T[2][2][1][0];
	A[2][5] = T[2][2][1][2];
	A[2][6] = T[2][2][2][1];
	A[2][7] = T[2][2][2][0];
	A[2][8] = T[2][2][0][2];
	
	A[3][0] = T[0][1][0][0];
	A[3][1] = T[0][1][1][1];
	A[3][2] = T[0][1][2][2];
	A[3][3] = T[0][1][0][1];
	A[3][4] = T[0][1][1][0];
	A[3][5] = T[0][1][1][2];
	A[3][6] = T[0][1][2][1];
	A[3][7] = T[0][1][2][0];
	A[3][8] = T[0][1][0][2];
	
	A[4][0] = T[1][0][0][0];
	A[4][1] = T[1][0][1][1];
	A[4][2] = T[1][0][2][2];
	A[4][3] = T[1][0][0][1];
	A[4][4] = T[1][0][1][0];
	A[4][5] = T[1][0][1][2];
	A[4][6] = T[1][0][2][1];
	A[4][7] = T[1][0][2][0];
	A[4][8] = T[1][0][0][2];
	
	A[5][0] = T[1][2][0][0];
	A[5][1] = T[1][2][1][1];
	A[5][2] = T[1][2][2][2];
	A[5][3] = T[1][2][0][1];
	A[5][4] = T[1][2][1][0];
	A[5][5] = T[1][2][1][2];
	A[5][6] = T[1][2][2][1];
	A[5][7] = T[1][2][2][0];
	A[5][8] = T[1][2][0][2];
	
	A[6][0] = T[2][1][0][0];
	A[6][1] = T[2][1][1][1];
	A[6][2] = T[2][1][2][2];
	A[6][3] = T[2][1][0][1];
	A[6][4] = T[2][1][1][0];
	A[6][5] = T[2][1][1][2];
	A[6][6] = T[2][1][2][1];
	A[6][7] = T[2][1][2][0];
	A[6][8] = T[2][1][0][2];
	
	A[7][0] = T[2][0][0][0];
	A[7][1] = T[2][0][1][1];
	A[7][2] = T[2][0][2][2];
	A[7][3] = T[2][0][0][1];
	A[7][4] = T[2][0][1][0];
	A[7][5] = T[2][0][1][2];
	A[7][6] = T[2][0][2][1];
	A[7][7] = T[2][0][2][0];
	A[7][8] = T[2][0][0][2];
	
	A[8][0] = T[0][2][0][0];
	A[8][1] = T[0][2][1][1];
	A[8][2] = T[0][2][2][2];
	A[8][3] = T[0][2][0][1];
	A[8][4] = T[0][2][1][0];
	A[8][5] = T[0][2][1][2];
	A[8][6] = T[0][2][2][1];
	A[8][7] = T[0][2][2][0];
	A[8][8] = T[0][2][0][2];
	
	
}

