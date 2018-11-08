#include "Matlib.h"
void zeroV(VectorD & A)
{
	fill(A.origin(), A.origin()+A.num_elements(), 0.0);
}

void zeroM2(MatrixD2 & A)
{
	fill(A.origin(), A.origin()+A.num_elements(), 0.0);
}

void zeroM3(MatrixD3 & A)
{
	fill(A.origin(), A.origin()+A.num_elements(), 0.0);
}

void zeroT2(Tensor2 & A)
{
	fill(A.origin(), A.origin()+A.num_elements(), 0.0);
}

void zeroT4(Tensor4 & A)
{
	fill(A.origin(), A.origin()+A.num_elements(), 0.0);
}

void onem(Tensor2 & A)
{
	zeroT2(A);
	
	for (int i = 0; i < 3; ++i ) A[i][i] = 1.0;
	
}

// Determinant of a 3x3 matrix
double matlib_determinant(MatrixD2 & A)
{
	double det;

	// 0 1 2 = [0][0] [0][1] [0][2]
	// 3 4 5 = [1][0] [1][1] [1][2]
	// 6 7 8 = [2][0] [2][1] [2][2]
	det = A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
         -A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
         +A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
   
	return det;
}

// Inversion of a 3x3 matrix, which also
// returns the determinant as a by-product.
// If the determinant is close to zero,
// no inversion is performed.
double matlib_inverse(MatrixD2 & A, MatrixD2 & Ainv)
{
	double det,detinv;

	det = matlib_determinant(A);
	if (fabs(det) < 1.e-20) return 0.e0;

	detinv = 1./det;
  
	// 0 1 2 = [0][0] [0][1] [0][2]
	// 3 4 5 = [1][0] [1][1] [1][2]
	// 6 7 8 = [2][0] [2][1] [2][2]

	Ainv[0][0] = detinv*( A[1][1]*A[2][2]-A[1][2]*A[2][1]);
	Ainv[0][1] = detinv*(-A[0][1]*A[2][2]+A[0][2]*A[2][1]);
	Ainv[0][2] = detinv*( A[0][1]*A[1][2]-A[0][2]*A[1][1]);
	Ainv[1][0] = detinv*(-A[1][0]*A[2][2]+A[1][2]*A[2][0]);
	Ainv[1][1] = detinv*( A[0][0]*A[2][2]-A[0][2]*A[2][0]);
	Ainv[1][2] = detinv*(-A[0][0]*A[1][2]+A[0][2]*A[1][0]);
	Ainv[2][0] = detinv*( A[1][0]*A[2][1]-A[1][1]*A[2][0]);
	Ainv[2][1] = detinv*(-A[0][0]*A[2][1]+A[0][1]*A[2][0]);
	Ainv[2][2] = detinv*( A[0][0]*A[1][1]-A[0][1]*A[1][0]);

	return det;
}

void matlib_transpose(MatrixD2 & A, MatrixD2 & AT)
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
