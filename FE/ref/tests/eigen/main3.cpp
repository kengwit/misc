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
	MatD2 E_TR(3,3);
	MatD2 epsE(3,3);
	MatD2 Dalg(6,6);
	MatD2 Sig(3,3);
	
	E_TR = 1.0,0.5,0.5,
		   0.5,1.0,0.5,
	       0.5,0.5,2.0;
	
	cout << "E_TR = " << endl << E_TR << endl;
	double tol = 1.e-6;
	double E,nu,sigmay;
	double Kmod,Gmod;
	double trSig,j2,f;
	EgVecD bm1; bm1.resize(9); bm1.setZero();
	
	EgVecD epsE6; epsE6.resize(6); 
	EgVecD epsE9; epsE9.resize(9); 
	EgVecD Sig6;  Sig6.resize(6);
	EgVecD Sig9;  Sig9.resize(9);
	EgVecD Sdev9; Sdev9.resize(9); 
	
	EgMatD2 De6; De6.resize(6,6); De6.setZero();
	EgMatD2 De9; De9.resize(9,9); De9.setZero();
	EgMatD2 Ce9; Ce9.resize(9,9); Ce9.setZero();
	EgMatD2 PmatT; PmatT.resize(6,9); PmatT.setZero(); // matrix to convert between 9- and 6-element stuff
	PmatT << 1.0,  0,  0,  0,  0,  0,  0,  0,  0,
		 	   0,1.0,  0,  0,  0,  0,  0,  0,  0,
			   0,  0,1.0,  0,  0,  0,  0,  0,  0,
			   0,  0,  0,0.5,0.5,  0,  0,  0,  0,
			   0,  0,  0,  0,  0,0.5,0.5,  0,  0,
			   0,  0,  0,  0,  0,  0,  0,0.5,0.5;
	
	cout << "PmatT = " << endl << PmatT << endl;
	
	// Obtain material properties
    E      = 1.0;
	nu     = 0.3;
	sigmay = 1.0;

	Gmod = E/(2.0*(1.0+nu));
	Kmod = E*nu/((1.0+nu)*(1.0-2.0*nu)); // this is really lambda (Lame's constant)
	
	// identity tensor size-9 vector form
	bm1(0) = 1.0;
	bm1(1) = 1.0;
	bm1(2) = 1.0;
	
	// convert elastic trial strain tensor to a size-9 vector		 
	epsE9(0) = E_TR(0,0); // xx  11
	epsE9(1) = E_TR(1,1); // yy  22
	epsE9(2) = E_TR(2,2); // zz  33
	epsE9(3) = E_TR(0,1); // xy  12
	epsE9(4) = E_TR(1,0); // yx  21
	epsE9(5) = E_TR(1,2); // yz  23
	epsE9(6) = E_TR(2,1); // zy  32
	epsE9(7) = E_TR(2,0); // zx  31
	epsE9(8) = E_TR(0,1); // zx  13
	
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
	
	// Trial Cauchy stress
	Sig9 = De9*epsE9;
	
	// Trial deviatoric stress vector
	trSig = (Sig9(0)+Sig9(1)+Sig9(2))/3.0;
	Sdev9 = Sig9-trSig*bm1;
	
	// Compute J2=1/2*s_ij*s_ij
	j2 = 0.5*Sdev9.dot(Sdev9);
	
	// yield function value
	f = sqrt(3.0*j2)/sigmay-1.0;

	
	// initialize total strain and Dalg to values 
	// corresponding to the elastic case 
	epsE = E_TR;
	
	Sig6 = PmatT*Sig9;
	

	Sig(0,0) = Sig6(0);
	Sig(1,1) = Sig6(1);
	Sig(2,2) = Sig6(2);
	Sig(0,1) = Sig6(3);
	Sig(1,0) = Sig6(3);
	Sig(1,2) = Sig6(4);
	Sig(2,1) = Sig6(4);
	Sig(2,0) = Sig6(5);
	Sig(0,2) = Sig6(5);
	
	cout << "Sig6 = " << endl << Sig6 << endl;
	cout << "Sig = " << endl << Sig << endl;
	
	De6 = PmatT*De9*PmatT.transpose();
	
	for ( int i = 0; i < 6; ++i ) {
		for ( int j = 0; j < 6; ++j ) {
			Dalg(i,j) = De6(i,j);
		}
	}
	
	cout << "De9 = " << endl << De9 << endl;
	cout << "De6 = " << endl << De6 << endl;
	cout << "Dalg = " << endl << Dalg << endl;
	EgVecD sdev; sdev.resize(9);
	
	EgVecD df; df.resize(9);
	double sdevnorm;
	sdev = Sig9;
	sdevnorm = sqrt(sdev.dot(sdev));
	cout << "sdevnorm = " << sdevnorm << endl;
	cout << "sdevnorm = " << sdev.norm() << endl;
	
	df = sdev/sdevnorm; //sqrt(3.0/2.0)/sigmay*sdev/sdevnorm;
   
	cout << "Sig9 = " << endl << Sig9 << endl;
	cout << "df = " << endl << df << endl;
	cout << "dfnorm = " << df.norm() << endl;
	
	EgMatD2 Id9; Id9.resize(9,9); Id9.setZero();
	EgMatD2 Idev; Idev.resize(9,9); Idev.setZero();
	EgMatD2 ddf; ddf.resize(9,9);
	for ( int i = 0; i < 9; ++i ) Id9(i,i) = 1.0;
	
	Idev = Id9 - bm1*bm1.transpose()/3.0;
	
	cout << "Idev = " << endl << Idev << endl;
	j2 = 1.0;
	ddf = sqrt(3)/(2.0*sigmay)*( Idev/sqrt(j2) - sdev*sdev.transpose()/(2*pow(sqrt(j2),3.0)) );
	cout << "ddf = " << endl << ddf << endl;
	
	EgVecD v1; v1.resize(4); 
	v1(0) = 2;
	v1(1) = 1;
	v1(2) = 1;
	v1(3) = 1.5;
	cout << v1.head(3).norm() << endl;
	cout << sqrt(6) << endl;
}
