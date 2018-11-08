

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cassert>
using namespace std;
#include "boost/multi_array.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
typedef Eigen::MatrixXd EgMatD2;
typedef Eigen::VectorXd EgVecD;

typedef boost::multi_array<double, 4> Tensor4;
typedef boost::multi_array<double, 2> Tensor2;
typedef boost::multi_array<double, 3> MatrixD3;
typedef boost::multi_array<double, 2> MatrixD2;
typedef boost::multi_array<double, 1> VectorD;


void zeroM2(MatrixD2 & A);
void zeroT2(Tensor2 & A);
void zeroT4(Tensor4 & A);
void onem(Tensor2 & A);
double matlib_determinant(MatrixD2 & A);
double matlib_inverse(MatrixD2 & A, MatrixD2 & Ainv);
void matlib_transpose(MatrixD2 & A, MatrixD2 & AT);
void matlib_matmult(Tensor2 & A, Tensor2 & B, Tensor2 & C);


int main()
{
	
	
	Tensor2 epsEtr(boost::extents[3][3]);
	
	double E,nu,sigmay;
	
	// Obtain material properties
    E      = 1.0;
    nu     = 0.3;
	sigmay = 0.01;
	
	
	// current elastic trial strain tensor
	zeroT2(epsEtr);
	epsEtr[0][0] = 0.1;
	epsEtr[1][1] = 0.1;
	epsEtr[0][1] = 0.05;
	epsEtr[1][0] = 0.05;
	
		
	// ===============================================================================================	
	double tol = 1.e-12;
	double Kmod,Gmod;
	double sdevnorm,trSig,j2,f;
	EgVecD bm1; bm1.resize(9); bm1.setZero();
	
	EgVecD epsE9(9); 
	EgVecD epsEtr9(9); 
	EgVecD Sig9(9);
	EgVecD sdev(9); 
	EgVecD nbar(9);
	 
	EgMatD2 Dalg6(6,6); Dalg6.setZero();
	EgMatD2 Dalg9(9,9); Dalg9.setZero();
	EgMatD2 De6(6,6);   De6.setZero();
	EgMatD2 De9(9,9);   De9.setZero();
	EgMatD2 Ce9(9,9);   Ce9.setZero();
	EgMatD2 Id9(9,9);   Id9.setZero();
	EgMatD2 Idev(9,9);  Idev.setZero();
	EgMatD2 Pmat(9,6);  Pmat.setZero(); // matrix to convert between 9- and 6-element stuff
	
	// identity tensor size-9 vector form
	bm1(0) = 1.0;
	bm1(1) = 1.0;
	bm1(2) = 1.0;
	
	for ( int i = 0; i < 9; ++i ) {
		Id9(i,i) = 1.0;
	}
	
	Idev = Id9 - bm1*bm1.transpose()/3.0;
	
	cout << "start Idev = " << endl << Idev << endl;
	Pmat << 1.0,  0,  0,  0,  0,  0,
		      0,1.0,  0,  0,  0,  0,
			  0,  0,1.0,  0,  0,  0,
			  0,  0,  0,0.5,  0,  0,
			  0,  0,  0,0.5,  0,  0,
			  0,  0,  0,  0,0.5,  0,
			  0,  0,  0,  0,0.5,  0,
			  0,  0,  0,  0,  0,0.5,
			  0,  0,  0,  0,  0,0.5;
			   
	
	
	// some material parameters
    Gmod = E/(2.0*(1.0+nu));     // shear modulus
	Kmod = E/(3.0*(1.0-2.0*nu)); // bulk modulus
	
	
	// convert elastic trial strain tensor to a size-9 vector		 
	epsEtr9(0) = epsEtr[0][0]; // xx  11
	epsEtr9(1) = epsEtr[1][1]; // yy  22
	epsEtr9(2) = epsEtr[2][2]; // zz  33
	epsEtr9(3) = epsEtr[0][1]; // xy  12
	epsEtr9(4) = epsEtr[1][0]; // yx  21
	epsEtr9(5) = epsEtr[1][2]; // yz  23
	epsEtr9(6) = epsEtr[2][1]; // zy  32
	epsEtr9(7) = epsEtr[2][0]; // zx  31
	epsEtr9(8) = epsEtr[0][2]; // zx  13
	
	// initialize elastic strain to trial elastic strain
	epsE9 = epsEtr9;
	
	// Isotropic elasticity 9x9 matrix
	double C11 = Kmod+4.0*Gmod/3.0;
	double C12 = Kmod-2.0*Gmod/3.0;
	double TWOG = 2.0*Gmod;
	De9 << C11, C12, C12,    0,    0,    0,    0,    0,    0,
	       C12, C11, C12,    0,    0,    0,    0,    0,    0,
	       C12, C12, C11,    0,    0,    0,    0,    0,    0,
	         0,   0,   0, TWOG,    0,    0,    0,    0,    0,
		     0,   0,   0,    0, TWOG,    0,    0,    0,    0,
		     0,   0,   0,    0,    0, TWOG,    0,    0,    0,
			 0,   0,   0,    0,    0,    0, TWOG,    0,    0,
			 0,   0,   0,    0,    0,    0,    0, TWOG,    0,
			 0,   0,   0,    0,    0,    0,    0,    0, TWOG; 
		
		
	// elastic compliance 9x9 matrix	
	Ce9 = De9.inverse();	
	
	//cout << "De9 = " << endl << De9 << endl;
	//cout << "Ce9 = " << endl << Ce9 << endl;
	
	
	// Trial Cauchy stress
	Sig9 = De9*epsEtr9;
	
	// Trial deviatoric stress vector
	trSig = (Sig9(0)+Sig9(1)+Sig9(2))/3.0;
	sdev = Sig9-trSig*bm1;
	
	// Compute J2=1/2*s_ij*s_ij
	sdevnorm = sdev.norm();
	j2 = 0.5*sdevnorm*sdevnorm;
	
	// yield function value
	f = sqrt(3.0*j2)/sigmay-1.0;
	//cout << "f = " << f << endl;
	//f = -1;
	if ( f > tol )
	{
		// ===============================
		// PLASTIC UPDATE
		// ===============================
		cout << "Plastic update\n";
		EgVecD resid; resid.resize(10); 
		EgVecD    df; df.resize(9);
		EgMatD2  ddf; ddf.resize(9,9);
		EgMatD2  dgam_ddf_De9; dgam_ddf_De9.resize(9,9);
		
		EgMatD2    A; A.resize(10,10);
		EgMatD2    Ainv; Ainv.resize(10,10);
		EgVecD    dx; dx.resize(10);
		
		int    maxit;
		double dgam;
		bool   converged_flag;
		
		converged_flag = false;
		maxit = 5;
		dgam  = 0.0;   // plastic multiplier
		
		resid.setZero(); // initialize residual vector
		resid(9)  = f;
		
		// =====================================
		// terms in the jacobian matrix
		// =====================================
		// notes:
		// d (j2) / dsigma = s
		// nvec x nvec = s/|s| x s/|s|, |s|^2 = 2*j2
		// d (|s|)^(-1) / dsigma
		// = d (2*j2)^(-1/2) / dsigma
		// = 0.5 * (2*j2)^(-3/2) * d (2*j2) / dsigma
		// = (2*j2)^(-3/2) * d (j2) / dsigma
		// = (2*j2)^(-3/2) * s
		nbar = sdev/sdevnorm;
		df   = sqrt(3.0/2.0)/sigmay*nbar;
		ddf  = sqrt(3.0)/(2.0*sigmay*sqrt(j2))*( Idev - nbar*nbar.transpose() ); 
   
		// newton iteration
		for ( int itnum = 1; itnum <= maxit; ++itnum )
		{
			cout << "\tnewton iter = " << itnum << " -------------------- " << endl;
			
			// form Jacobian
			A.setZero();
			A.topLeftCorner(9,9) = Id9+dgam*ddf*De9;
			A.block(9,0,1,9) = df.transpose()*De9;
			A.block(0,9,9,1) = df;
			
			//cout << "Jacobian A = " << endl << A << endl;
			//cin.get();
			
			Ainv = A.inverse();
			
			dx = Ainv * resid;
			
			//cout << "dx = " << endl << dx << endl;
			//cin.get();
			
			for ( int i = 0; i <=8; ++i ) {
				epsE9(i) -= dx(i); // elems 1 through 9
			}			
			dgam  -= dx(9);      // elem 10
			Sig9 = De9*epsE9;
			
			cout << "Sig9 = " << endl << Sig9 << endl;
			
			// deviatoric stress vector
			trSig = (Sig9(0)+Sig9(1)+Sig9(2))/3.0;
			sdev = Sig9-trSig*bm1;
	
			// Compute J2=1/2*s_ij*s_ij
			sdevnorm = sdev.norm();
			j2 = 0.5*sdevnorm*sdevnorm;
			
			cout << "j2 = " << endl << j2 << endl;
	
			// recompute terms in jacobian matrix for next iteration
			nbar = sdev/sdevnorm;
			df   = sqrt(3.0/2.0)/sigmay*nbar;
			cout << "nbar*nbar.transpose() = " << endl << nbar*nbar.transpose() << endl;
			cout << "sqrt(3.0)/(2.0*sigmay*sqrt(j2)) = " << sqrt(3.0)/(2.0*sigmay*sqrt(j2)) << endl;
			cout << "Idev = " << endl << Idev << endl;
			ddf  = sqrt(3.0)/(2.0*sigmay*sqrt(j2))*( Idev - nbar*nbar.transpose() ); 
   			
   			cout << "df = " << endl << df << endl;
   			cout << "ddf = " << endl << ddf << endl;
   			
			// recompute residual vector
			resid.head(9) = epsE9 - epsEtr9 + dgam*df;
			resid(9) = sqrt(3.0*j2)/sigmay-1.0;
			
			
			if ( ( resid.head(9).norm() < tol ) && ( abs(resid(9)) < tol ) )
			{
				converged_flag = true;
				break;
			}
			
		} // end return mapping 
		
		if ( converged_flag == true ) {
			
			
			cout << "df = " << endl << df << endl;
   			cout << "ddf = " << endl << ddf << endl;
			cout << "dgam = " << endl << dgam << endl;
			
			
			// consistent tangent operator
			A.setZero();
			A.topLeftCorner(9,9) = Ce9 + dgam*ddf;
			A.block(9,0,1,9) = df.transpose();
			A.block(0,9,9,1) = df;
							
			EgMatD2 B = A.inverse();
			cout << "Ainv = " << endl << Ainv << endl;
			
			Dalg9 = B.topLeftCorner(9,9); 
			Dalg6 = Pmat.transpose()*Dalg9*Pmat;
			
		} else {
			cout << "no convergence in return map after " << maxit << " iterations " << endl;
			exit(1);	
		}
		
	} else {
		
		// ===============================
		// ELASTIC UPDATE
		// ===============================
		cout << "Elastic update\n";
		De6 = Pmat.transpose()*De9*Pmat;	
		Dalg6 = De6;
		

	} // end check yield	
	
	cout << "Dalg6 = " << endl << Dalg6 << endl;
	
	return 0;
	
}

// -------------------------------------------------
// CONVENIENCE FUNCTIONS
// -------------------------------------------------
void zeroM2(MatrixD2 & A)
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
	for (int i = 0; i < 3; ++i )
		for (int j = 0; j < 3; ++j )
			A[i][j] = 0.0;
	
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

/*  

To link the object file with ABAQUS, the command link_sl of the ABAQUS environmental file abaqus_v6.env must be modified so that the C++ libraries are correctly loaded by the Fortran compiler, which is used for linking. The correct link_sl command is the following:

link_sl='LINK', '/nologo', '/INCREMENTAL:NO', '/subsystem:console', '/machine:X64','/DEFAULTLIB:DFORMD.LIB', '/NODEFAULTLIB:LIBC.LIB', '/DEFAULTLIB:OLDNAMES.LIB', '/DEFAULTLIB:MSVCRT.LIB', '/NODEFAULTLIB:LIBCMT.LIB', '/NODEFAULTLIB:LIBCPMT.LIB', '/DEFAULTLIB:LIBCMT.LIB', '/DEFAULTLIB:kernel32.lib', '/DEFAULTLIB:user32.lib', '/DEFAULTLIB:advapi32.lib', '/FIXED:NO', '/dll','/NOENTRY', '/def:%E', '/out:%U', '%F', '%A', '%B', 'oldnames.lib', 'user32.lib', 'ws2_32.lib', 'netapi32.lib', 'advapi32.lib' 

link_sl='LINK', '/nologo', '/INCREMENTAL:NO', '/subsystem:console', '/machine:I386','/DEFAULTLIB:DFORMD.LIB', '/NODEFAULTLIB:LIBC.LIB', '/DEFAULTLIB:OLDNAMES.LIB', '/DEFAULTLIB:MSVCRT.LIB', '/NODEFAULTLIB:LIBCMT.LIB', '/NODEFAULTLIB:LIBCPMT.LIB', '/DEFAULTLIB:LIBCMT.LIB', '/DEFAULTLIB:kernel32.lib', '/DEFAULTLIB:user32.lib', '/DEFAULTLIB:advapi32.lib', '/FIXED:NO', '/dll','/NOENTRY', '/def:%E', '/out:%U', '%F', '%A', '%B', 'oldnames.lib', 'user32.lib', 'ws2_32.lib', 'netapi32.lib', 'advapi32.lib'

Then, we may start the Abaqus job with the DOS command line

abaqus job=inputfile user=umat_vec.obj

compile_cpp=['cl', '/c', '/nologo', '/W0', '/MD', 
             '/TP', '/EHsc', '/DNDEBUG', '/DWIN32', '/DTP_IP', '/D_CONSOLE', 
             '/DNTI', '/DFLT_LIC', '/DOL_DOC', '/D__LIB__', '/DHKS_NT',
             '/DABQ_NTI_NET', '/DFAR=', '/D_WINDOWS', '/DABQ_WIN86_64', 
             '/O1', '/I%I']

compile_fortran=['ifort', '/c','/DABQ_WIN86_64',
                 '/iface:cref', '/recursive', '/Qauto-scalar', '/nologo',
                 '/heap-arrays:1', '/Od', '/include:%I']

link_sl=['LINK', '/nologo', '/INCREMENTAL:NO', '/subsystem:console', '/machine:AMD64', '/NODEFAULTLIB:LIBC.LIB', '/NODEFAULTLIB:LIBCMT.LIB', 
         '/DEFAULTLIB:OLDNAMES.LIB', '/DEFAULTLIB:LIBIFCOREMD.LIB', '/DEFAULTLIB:LIBIFPORTMD', '/DEFAULTLIB:LIBMMD.LIB',
         '/DEFAULTLIB:MSVCRT.LIB', '/DEFAULTLIB:kernel32.lib', '/DEFAULTLIB:user32.lib', '/DEFAULTLIB:advapi32.lib',
         '/FIXED:NO', '/dll', '/def:%E', '/out:%U', '%F', '%A', '%L', '%B', 'oldnames.lib', 'user32.lib', 'ws2_32.lib',
         'netapi32.lib', 'advapi32.lib']
#, '&&', 'mt', '/manifest', '%U.manifest', '/outputresource:%U', '&&', 'del', '%U.manifest']

link_exe=['LINK', '/nologo', '/INCREMENTAL:NO', '/subsystem:console', '/machine:AMD64', '/STACK:20000000',
          '/NODEFAULTLIB:LIBC.LIB', '/NODEFAULTLIB:LIBCMT.LIB', '/DEFAULTLIB:OLDNAMES.LIB', '/DEFAULTLIB:LIBIFCOREMD.LIB',
          '/DEFAULTLIB:LIBIFPORTMD', '/DEFAULTLIB:LIBMMD.LIB', '/DEFAULTLIB:MSVCRT.LIB', '/DEFAULTLIB:kernel32.lib',
          '/DEFAULTLIB:user32.lib', '/DEFAULTLIB:advapi32.lib', '/FIXED:NO', '/LARGEADDRESSAWARE', '/out:%J', '%F', '%M',
          '%L', '%B', '%O', 'oldnames.lib', 'user32.lib', 'ws2_32.lib', 'netapi32.lib', 'advapi32.lib']
#, '&&', 'mt', '/manifest', '%J.manifest', '/outputresource:%J', '&&', 'del', '%J.manifest']

# Link command to be used for MAKE w/o fortran compiler.
# remove the pound signs in order to remove the comments and have the file take effect.
#
#link_exe=['LINK', '/nologo', 'INCREMENTAL:NO', '/subsystem:console', '/machine:AMD64', '/NODEFAULTLIB:LIBC.LIB', '/NODEFAULTLIB:LIBCMT.LIB',
#          '/DEFAULTLIB:OLDNAMES.LIB', '/DEFAULTLIB:MSVCRT.LIB', '/DEFAULTLIB:kernel32.lib', 'DEFAULTLIB:user32.lib', '/DEFAULTLIB:advapi32.lib',
#          '/FIXED:NO', '/LARGEADDRESSAWARE', '/out:%J', '%F', '%M', '%L', '%B', '%O', 'oldnames.lib', 'user32.lib', 'ws2_32.lib',
#          'netapi32.lib', 'advapi32.lib', '&&', 'mt', '/manifest', '%J.manifest', '/outputresource:%J', '&&', 'del', '%J.manifest']

usub_lib_dir='C:/Users/KWL/Desktop/Research/abaqus_umats/kengwit/'
 
*/

/*
This array

    \begin{bmatrix} 11 & 12 & 13 \\ 
	                21 & 22 & 23 \end{bmatrix}

Would be stored as follows in the two orders:
Column-Major Order
e.g. Fortran Address 	Value
0 	11
1 	21
2 	12
3 	22
4 	13
5 	23
*/

