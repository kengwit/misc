#ifndef NEOHOOKEAN_H
#define NEOHOOKEAN_H

#include "Material.h"

typedef boost::multi_array<double, 2> Tensor2;
typedef boost::multi_array<double, 4> Tensor4;


class NeoHookean : public Material
{
public:
	NeoHookean() { 
		//cout << "NeoHookean construct" << endl;
		mat_type  = "NeoHookean";
		//
		// number of properties = 2
		// lambda 
		// mu     
		//
		nprops  = 2; 
		nsvars  = 9; // we will store Fe_old as state variables at each gauss point of element
		
		allocate_data();		
		
	} 
    ~NeoHookean() { 
		//cout << "LinearElastic destroy" << endl; 		
		clear_data();
	} 
	
	void initialize_state_variables(ViewD1 & svars);
    void constitutive_update(unsigned int el, unsigned int gp, unsigned int sec_type, VectorD & props, ViewD2 & eps_t, ViewD2 & eps_tau, ViewD2 & Dmat, ViewD1 & Svec, ViewD1 & svars, bool tangent_flag); 
};


// Material utility functions
template <typename T>
void zero(T & A)
{
	fill(A.origin(), A.origin()+A.num_elements(), 0.0);
}

template <typename T>
void onem(T & A)
{
	zero(A);
	
	for (int i = 0; i < 3; ++i ) A[i][i] = 1.0;
	
}

// Determinant of a 3x3 matrix
template <typename T>
double matlib_determinant(T & A)
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
template <typename T1, typename T2>
double matlib_inverse(T1 & A, T2 & Ainv)
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

template <typename T1, typename T2>
void matlib_transpose(T1 & A, T2 & AT)
{
	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
			AT[i][j] = A[j][i];
}


template <typename T1, typename T2, typename T3>
void matlib_matmult(T1 & A, T2 & B, T3 & C)
{
	zero(C);
	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
			for ( int k = 0; k < 3; ++k )
				C[i][j] += A[i][k]*B[k][j];
}



#endif