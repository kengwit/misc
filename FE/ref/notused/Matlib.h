#ifndef MATLIB_H
#define MATLIB_H

#include <cmath>
using namespace std;

#include "Globals.h"

	
void zeroV(VectorD & A);
void zeroM2(MatrixD2 & A);
void zeroM3(MatrixD3 & A);
void zeroT2(Tensor2 & A);
void zeroT4(Tensor4 & A);
double matlib_determinant(MatrixD2 & A);
double matlib_inverse(MatrixD2 & A, MatrixD2 & Ainv);
void matlib_transpose(MatrixD2 & A, MatrixD2 & AT);
void matlib_matmult(Tensor2 & A, Tensor2 & B, Tensor2 & C);

	
	
#endif