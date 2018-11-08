/*!************************************************************************
! User element (UEL) for XXX material -- plane-strain
!************************************************************************
! Element details:
!************************************************************************
!
! Solution variables (or nodal variables) are the displacements (DOFs 1-2)
!
! Material behavior is XXX material 
! 
! This subroutine is for a two-dimensional 4-node isoparametric
!  quadrilateral element as shown below with 4pt integration with F-bar method 
!
!
!              A eta (=xi_2)
!  4-node      |
!   quad       | Face 3
!        4-----------3
!        |     |     |
!        |     |     |
! Face 4 |     ------|---> xi (=xi_1)
!        |           |
!        |           |  Face 2
!        1-----------2
!            Face 1
!
!
!************************************************************************
! Usage:
!************************************************************************
!
! User element statement in the input file:
!  *User Element,Nodes=4,Type=U1,Iproperties=1,Properties=3,Coordinates=2,Variables=9,Unsymm
!  1,2
!
! Note: 36 state variables are used in this element to store previously 
! converged elastic deformation gradient (0.5*log(Be))
!
! NODES = 4
! TYPE=U1 (plane strain)
! IPROPERTIES = 9 (number of state vars per integration point)
! PROPERTIES = E,nu,sigmay (3 properties for elastic perfectly plastic)
! COORDINATES = 2 at each node
! VARIABLES = 4*9 = 36 (Store Fe_tau at each integration point)
! UNSYMM = use unsymmetric solver for F-bar method
! All dofs are displacements 1,2
!
The C++ umat must be compiled externally using a compiler supported by ABAQUS, producing an object file .obj. The compiler command is

cl /EHsc /c /Fo uel4nodeFbarEPBoostStoreFeUntested.cpp /I"C:\Users\KWL\Desktop\Research\libraries\boost_1_57_0" /I"C:\Users\KWL\Desktop\Research\libraries\eigen-3.2.4" 
cl /O2 /DNDEBUG /EHsc /c /Fo uel4nodeFbarEPBoostStoreFeUntested.cpp /I"C:\Users\KWL\Desktop\Research\libraries\boost_1_57_0" /I"C:\Users\KWL\Desktop\Research\libraries\eigen-3.2.4" 

Then do:
abaqus job=XXX user=uel4nodeFbarXXXXX.obj interactive

!*************************************************************************/


//#include <aba_for_c.h> since 6.13?
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cassert>
#include <limits>
using namespace std;

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
typedef Eigen::MatrixXd EgMatD2;
typedef Eigen::VectorXd EgVecD;

//#define BOOST_DISABLE_ASSERTS 
#include "boost/multi_array.hpp"
typedef boost::multi_array<double, 4> Tensor4;
typedef boost::multi_array<double, 2> Tensor2;
typedef boost::multi_array<double, 3> MatrixD3;
typedef boost::multi_array<double, 2> MatrixD2;
typedef boost::multi_array<double, 1> VectorD;

// column-major for fortran array access
// h = row-height
#define FORT_ACCESS(A,r,c,h) ((A)[(r) + (c)*(h)])

struct eigen_struct { 
	int index;
	double value;
	double vector[3];
	bool operator<(eigen_struct const &other) const { 
		return value < other.value;
	}
};
void solve_eigenvalue_3x3(Tensor2 & AMat, VectorD & APr, MatrixD2 & Nvec);

void zeroV(VectorD & A);
void zeroM2(MatrixD2 & A);
void zeroM3(MatrixD3 & A);
void zeroT2(Tensor2 & A);
void zeroT4(Tensor4 & A);
void onem(Tensor2 & A);
double matlib_determinant(Tensor2 & A);
double matlib_inverse(Tensor2 & A, Tensor2 & Ainv);
void matlib_transpose(Tensor2 & A, Tensor2 & AT);
void matlib_matmult(Tensor2 & A, Tensor2 & B, Tensor2 & C);
void convertTensor4ToMat9(Tensor4 & T, MatrixD2 & A);
void convertMat9ToTensor4(MatrixD2 & A, Tensor4 & T);

void UPE4(double * RHS, double * AMATRX, double * SVARS, double * ENERGY, int * NDOFEL, int * NRHS, int * NSVARS, double * PROPS, int * NPROPS, 
		  double * COORDS, int * MCRD, int * NNODE, double * U, double * DU, double * V, double * A, int * JTYPE, double * TIME, double * DTIME,
          int * KSTEP, int * KINC, int * JELEM, double * PARAMS, int * NDLOAD, int * JDLTYP, double * ADLMAG, double * PREDEF,
          int * NPREDF, int * LFLAGS,int * MLVARX,double * DDLMAG,int * MDLOAD,double * PNEWDT,int * JPROPS,int * NJPROP, double * PERIOD,
		  int nquad, int ndm);
		  
void compute_weight_upe4(int nintpts, VectorD & w);
void compute_shape_upe4(int ndm, int nen, int nintpts, MatrixD2 & shp, MatrixD3 & dshp, VectorD & detj, MatrixD2 & x, bool shp_flag);

void constitutive_update(double * PROPS, Tensor2 & epsE_t, Tensor2 & F_t, Tensor2 & F_tau,
						 Tensor2 & epsE_tau, Tensor2 & T_tau, Tensor4 & SpTanMod);

void VMconst(double * PROPS, Tensor2 & epsEtr,
			 MatrixD2 & Dalg, Tensor2 & sig_tau, Tensor2 & epsE_tau);

void parDerGen(Tensor2 & X, VectorD & XPr, MatrixD2 & Nvec_X, Tensor4 & LTensor);
			 
// note function umat with underscore "uel_" does not work
extern "C" void uel(double * RHS, double * AMATRX, double * SVARS, double * ENERGY, int * NDOFEL, int * NRHS, int * NSVARS, double * PROPS, int * NPROPS, 
                    double * COORDS, int * MCRD, int * NNODE, double * U, double * DU, double * V, double * A, int * JTYPE, double * TIME, double * DTIME,
					int * KSTEP, int * KINC, int * JELEM, double * PARAMS, int * NDLOAD, int * JDLTYP, double * ADLMAG, double * PREDEF,
					int * NPREDF, int * LFLAGS,int * MLVARX,double * DDLMAG,int * MDLOAD,double * PNEWDT,int * JPROPS,int * NJPROP, double * PERIOD)
{
	
	// Check the procedure type, this should be a 
	// static analysis, which is either 1 or 2
    if((LFLAGS[0]==1)||(LFLAGS[0]==2)) {	
        // correct procedure specified
    } else {
        cout << "Abaqus does not have the right procedure" << endl;
        cout << "go back and check the procedure type" << endl;
        cout << "LFLAGS[0]=" << LFLAGS[0] << endl;
        exit(1);
    }


    // Make sure Abaqus knows you are doing a large
    // deformation problem
    if(LFLAGS[1]==0) {	
        // LFLAGS[1]=0 -> small disp.
        // LFLAGS[1]=1 -> large disp.
        cout << "Abaqus thinks you are doing" << endl;
        cout << "a small displacement analysis" << endl;
        cout << "go in and set nlgeom=yes" << endl;
        exit(1);
    }

	
    // Check to see if you are doing a general
    //  step or a linear purturbation step
    if(LFLAGS[3]==1) {	
        // LFLAGS[3]=0 -> general step
        // LFLAGS[3]=1 -> linear purturbation step
        cout << "Abaqus thinks you are doing" << endl;
        cout << "a linear purturbation step" << endl;
        exit(1);
	}		 
 
    // Do nothing if a ``dummy'' step
    if( *DTIME == 0.0 ) return;
	
	//----------------------------------------------------------------
    // Call the particular element to perform the analysis
    //----------------------------------------------------------------
    if ( *JTYPE == 1 ) { // *JTYPE = the number you used to label the user element e.g. U1 for plane strain
		int nquad =  4;
		int ndm   =  2;
		
		UPE4(RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS, PROPS, NPROPS, 
		     COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, TIME, DTIME,
			 KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG, PREDEF,
			 NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT, JPROPS, NJPROP, PERIOD,
			 nquad, ndm);
	
	} else {
		cout << "unrecognized element\n";
	}
      
}

void compute_weight_upe4(int nintpts, VectorD & w)
{
	assert ( nintpts == 1 || nintpts == 4 );
	
	if ( nintpts == 1 ) {
		
		w[0] = 2*2;
		
	} else if ( nintpts == 4 ) {
	
		w[0] = 1.0*1.0;
        w[1] = 1.0*1.0;
		w[2] = 1.0*1.0; 
	    w[3] = 1.0*1.0;
		
	} else {
		cout << "compute_weight_upe4: shouldn't be here\n";	
	}
	
	
}

void compute_shape_upe4(int ndm, int nen, int nintpts, MatrixD2 & shp, MatrixD3 & dshp, VectorD & detj, MatrixD2 & x, bool shp_flag)
{
	
	assert ( nintpts == 1 || nintpts == 4 );
	
	MatrixD2    s(boost::extents[nintpts][2]);
	MatrixD2 dnds(boost::extents[nen][ndm]); // temp
	MatrixD2 xjac(boost::extents[ndm][ndm]); // temp
	MatrixD2 dxds(boost::extents[ndm][ndm]); // temp
	double xsj;           		             // temp
	
	double forth = 0.25;
	
	if ( nintpts == 1 ) {
		
		s[0][0] = 0.0;
		s[0][1] = 0.0;	
	
	} else if ( nintpts == 4 )  {
		
		//long double pt = sqrt(static_cast<long double>(1.0)/static_cast<long double>(3.0));
		double pt = sqrt(1.0/3.0);
		//cout << "pt = ";
		//cout << std::fixed;
		//cout.precision(51);
		//cout << pt << endl;
			
		s[0][0] = -pt;
		s[0][1] = -pt;
		
		s[1][0] =  pt;
		s[1][1] = -pt;
		
		s[2][0] =  pt;
		s[2][1] =  pt;
		
		s[3][0] = -pt;
		s[3][1] =  pt;
		
	} else {
	
		cout << "compute_shape_upe4: shouldn't be here\n";
	
	}

	if ( shp_flag == true ) {
		zeroM2(shp);
	}
	
	zeroM3(dshp);
	zeroV(detj);
	
	for (int lquad = 0; lquad < nintpts; lquad++)
	{
		//cout << "gauss coord = (" << s[lquad][0] << "," << s[lquad][1] << endl;
		// ===================================
		// fill in shape functions N1,N2,N3,N4
		// ===================================
		if ( shp_flag == true ) {
			shp[lquad][0] = forth*(1.0-s[lquad][0])*(1.0-s[lquad][1]); // N1 
			shp[lquad][1] = forth*(1.0+s[lquad][0])*(1.0-s[lquad][1]); // N2 
			shp[lquad][2] = forth*(1.0+s[lquad][0])*(1.0+s[lquad][1]); // N3 
			shp[lquad][3] = forth*(1.0-s[lquad][0])*(1.0+s[lquad][1]); // N4 		
		}
		
		// ===============================
		// Compute shape function derivatives wrt natural coordinates 
		//  - this is an intermediate stage
		// ===============================
		// shape function derivatives (d/ds)
		// For e.g.:
		// dN/ds = [  dN1/dr   dN1/ds;
		//                 ....
		//          dNnen/dr dNnen/ds]; nen*ndm
		// 
		dnds[0][0] = -forth*(1.0-s[lquad][1]); // dN1/dr
		dnds[0][1] = -forth*(1.0-s[lquad][0]); // dN1/ds

		dnds[1][0] =  forth*(1.0-s[lquad][1]); // dN2/dr
		dnds[1][1] = -forth*(1.0+s[lquad][0]); // dN2/ds

		dnds[2][0] =  forth*(1.0+s[lquad][1]); // dN3/dr
		dnds[2][1] =  forth*(1.0+s[lquad][0]); // dN3/ds
		
		dnds[3][0] = -forth*(1.0+s[lquad][1]); // dN4/dr
		dnds[3][1] =  forth*(1.0-s[lquad][0]); // dN4/ds
		
		// ==============================
		// Compute Jacobian matrix
		// ==============================
		// 3D dxds = [dx/dr dy/dr dz/dr; = [0 1 2 
		//            dx/ds dy/ds dz/ds;    3 4 5
		//			  dx/dt dy/dt dz/dt]    6 7 8];
		//
		// 2D dxds = [dx/dr dy/dr;  = [0 1
		//            dx/ds dy/ds];    2 3];
		//
		// 1D dxds = [dx/ds] 
		// 
		// (0,0) dx/dr = dN1/dr*x1 + ... + dN3/dr*x3 = 1*0+0*1-1*0 = 0
		// (0,1) dy/dr = dN1/dr*y1 + ... + dN3/dr*y3 
		// (1,0) dx/ds = dN1/ds*x1 + ... + dN3/ds*x3 
		// (1,1) dy/ds = dN1/ds*y1 + ... + dN3/ds*y3 	
		
		zeroM2(dxds);
		for ( int i = 0; i < ndm; ++i ) {
			for ( int j = 0; j < ndm; ++j ) {
				for ( int k = 0; k < nen; ++k ) {
					dxds[i][j] += dnds[k][i]*x[k][j];
				}
			}
		}
		
		// ==============================
		// Compute determinant of Jacobian matrix
		// ==============================
		xsj = dxds[0][0]*dxds[1][1]-dxds[1][0]*dxds[0][1];
		//cout << "xsj = " << xsj << endl;
		assert( xsj > 1.e-20 );
		
		// ==============================
		// Compute INVERSE of Jacobian matrix
		// ==============================
		// note: inverse of Jac is NOT just the Jac with terms inversed!
		// 3D xjac = [dr/dx ds/dx dt/dx; = [0 1 2
		//            dr/dy ds/dy dt/dy     3 4 5
		//            dr/dz ds/dz dt/dz];   6 7 8];
		//
		// 2D xjac = [dr/dx ds/dx;  = [0 1
		//            dr/dy ds/dy];    2 3];
		//
		// 1D xjac = [ds/dx];
		//		
		// below is xjac = cofactor(dxds)^Transpose / det(dxds)
		xjac[0][0] = dxds[1][1]/xsj; // dr/dx
		xjac[0][1] =-dxds[0][1]/xsj; // ds/dx
		xjac[1][0] =-dxds[1][0]/xsj; // dr/dy
		xjac[1][1] = dxds[0][0]/xsj; // ds/dy
		
		// ==============================
		// store determinant of Jacobian 
		// matrix times integration weight
		// ==============================
		detj[lquad] = xsj; // xsj=ratio of real space area to reference area (4 units) e.g. quarter for unit square in real space
    	
		// ==============================
		// compute and store shape function 
		// derivatives wrt to real coordinates
		// ==============================
		// just think of Jac^inv * {dNi/dr}
		//                         {dNi/ds} 	
		for ( int a = 0; a < nen; ++a )
			for ( int i = 0; i < ndm; ++i )
				for ( int j = 0; j < ndm; ++j )
					dshp[lquad][a][i] += xjac[i][j]*dnds[a][j];
	}

}

void UPE4(double * RHS, double * AMATRX, double * SVARS, double * ENERGY, int * NDOFEL, int * NRHS, int * NSVARS, double * PROPS, int * NPROPS, 
		  double * COORDS, int * MCRD, int * NNODE, double * U, double * DU, double * V, double * A, int * JTYPE, double * TIME, double * DTIME,
          int * KSTEP, int * KINC, int * JELEM, double * PARAMS, int * NDLOAD, int * JDLTYP, double * ADLMAG, double * PREDEF,
          int * NPREDF, int * LFLAGS,int * MLVARX,double * DDLMAG, int * MDLOAD, double * PNEWDT, int * JPROPS,int * NJPROP, double * PERIOD, 
		  int nquad, int ndm)
{
	int counter;  		  
	int ndf   = ndm;
	int nsdm  = ndm*ndf;
	int nen   = *NNODE;
	int nlSdv = JPROPS[0];  // number of state variables
	
	Tensor2 Fbar_t(boost::extents[3][3]);   double detFbar_t;   // previous converged Fbar
	Tensor2 Fbar_tau(boost::extents[3][3]); double detFbar_tau; // current Fbar
	
	Tensor2 F_t(boost::extents[3][3]);    double detF_t;       // previous converged def gradient
	Tensor2 Fe_t(boost::extents[3][3]);                        // previous converged ELASTIC def gradient
	Tensor2 F_tau(boost::extents[3][3]);  double detF_tau;     // current def gradient
	Tensor2 Fe_tau(boost::extents[3][3]);                      // current ELASTIC def gradient
	Tensor2 Fc_tau(boost::extents[3][3]); double detFc_tau;    // current def gradient at centroid 'c' 
	Tensor2 Fc_t(boost::extents[3][3]);   double detFc_t;      // previous converged def gradient at centroid 'c' 
	Tensor2 T_tau(boost::extents[3][3]);
	Tensor4 SpTanMod(boost::extents[3][3][3][3]); 
	
	
	// ------------------------------------------------------
	// shape functions & derivatives at centroid
	// ------------------------------------------------------
	MatrixD2 shp0(boost::extents[1][nen]);         // only 1 gauss point (at centroid), not used for plane strain (used in G0mat for axisymm)
	
	MatrixD3 dshp0(boost::extents[1][nen][ndm]);   // ref config, used to calculate deformation gradient at centroid
	VectorD  detj0(boost::extents[1]);             // ref config, not used
	
	MatrixD3 dshpC0(boost::extents[1][nen][ndm]);  // current config, used to calculate G0mat
	VectorD  detjC0(boost::extents[1]);            // current config, not used
	
	// ------------------------------------------------------
	// shape functions & derivatives at normal gauss points
	// ------------------------------------------------------
	MatrixD2 shp(boost::extents[nquad][nen]);       // shape functions
	
	MatrixD3 dshp(boost::extents[nquad][nen][ndm]); // ref config, used to calculate deformation gradient
	VectorD  detj(boost::extents[nquad]);           // ref config, not used
		
	MatrixD3 dshpC(boost::extents[nquad][nen][ndm]); // current config, used for Bmat, Gmat
	VectorD  detjC(boost::extents[nquad]);           // current config, used for Bmat, Gmat
	  	
	// ---------
	// weights 
	// ---------  
	VectorD  w(boost::extents[nquad]);  // weights
	
	// ------------------------------------------------------
	// coordinates & displacements
	// ------------------------------------------------------
	MatrixD2 u(boost::extents[nen][ndm]);
	MatrixD2 du(boost::extents[nen][ndm]);
	MatrixD2 uOld(boost::extents[nen][ndm]);
	MatrixD2 coords(boost::extents[nen][ndm]);
	MatrixD2 coordsC(boost::extents[nen][ndm]);
	
	// ------------------------------------------------------
	// Voigt style quantities
	// ------------------------------------------------------
	MatrixD2 Kuu(boost::extents[nen*ndm][nen*ndm]);
	MatrixD2 Ru(boost::extents[nen][ndm]);
	VectorD  body(boost::extents[ndm]);
	
	VectorD  Smat(boost::extents[3]);
	MatrixD2 Bmat(boost::extents[3][ndm*nen]);
	MatrixD2 Gmat(boost::extents[4][ndm*nen]);
	MatrixD2 G0mat(boost::extents[4][ndm*nen]);
	MatrixD2 Qmat(boost::extents[4][4]);
	MatrixD2 Amat(boost::extents[4][4]);
	
	
    // Initialize the residual and tangent matrices to zero.
    zeroM2(Kuu);
    zeroM2(Ru);
    *ENERGY=0.0; 
	
	// Body forces
	zeroV(body);
	
    // Obtain nodal displacements 
	counter = 0;
	for ( int a = 0; a < nen; ++a ) {
		for ( int j = 0; j < ndm; ++j ) {
			// u[a][j]   =  U[j+a*ndm];
            //du[a][j]   = DU[j+a*ndm];
			 u[a][j]   =  U[counter];
            du[a][j]   = DU[counter];
			
			uOld[a][j] = u[a][j]-du[a][j];
			//cout << "U[" << counter << "] = " << U[counter] << endl;
			//cout << "DU[" << counter << "] = " << DU[counter] << endl;
			counter++;			
		}
	}
	
	// COORDS is column major
	// 0 2 4
	// 1 3 5
	counter = 0;
	for ( int a = 0; a < nen; ++a ) {
		for ( int j = 0; j < ndm; ++j ) {
			 //coords[i][j] = COORDS[FORT_ACCESS(COORDS,j,a,ndm);
			
			 coords[a][j] = COORDS[j+a*ndm];
			coordsC[a][j] = coords[a][j] + u[a][j];		
			//coordsC[a][j] = COORDS[j+a*ndm]+U[j+a*ndm];		
			
			//cout << "coords[" << i << "][" << j << "] = " << coords[i][j] 
			//	 << ", coordsC[" << i << "][" << j << "] = " << coordsC[i][j] << endl; 
			
		}
	}
	
	//cout << "area = " << coordsC[2][0]*coordsC[2][1] << endl;
	// ----------------------------------------------------------------
	// compute shape functions and derivatives
    if(nen == 4) {
		
		compute_shape_upe4(ndm, nen,     1, shp0,  dshp0,  detj0,  coords, true ); // at centroid (0), reference, for F-bar method
		compute_shape_upe4(ndm, nen,     1, shp0, dshpC0, detjC0, coordsC, false); // at centroid (0), current (C), for F-bar method
        compute_shape_upe4(ndm, nen, nquad,  shp,   dshp,   detj,  coords, true ); // at all gauss points, reference
        compute_shape_upe4(ndm, nen, nquad,  shp,  dshpC,  detjC, coordsC, false); // at all gauss points, current (C)
        compute_weight_upe4(nquad, w);
		
    } else {
		cout << "Incorrect number of nodes: nNode != 4" << endl;
        exit(1);
	}
   
	// ==================================================================
	// Calculate the deformation gradient at the element centriod
    // at the the begining and end of the increment for use in 
    // the `F-bar' method. `Tau' represents the end of the increment
    // and `t' the previous increment.
    onem(Fc_tau); // current
    onem(Fc_t);   // previous
	//zeroM2(Fc_tau);
	//zeroM2(Fc_t);
	
	for ( int i=0; i < ndm; ++i ) {
		for ( int j=0; j < ndm; ++j ) {
			for ( int a=0; a < nen; ++a ) {
				Fc_tau[i][j] +=    u[a][i]*dshp0[0][a][j];
				Fc_t[i][j]   += uOld[a][i]*dshp0[0][a][j];
			}
		}
	}
	
	
	// modify for plane-strain
	Fc_tau[2][2] = 1.0;
	Fc_t[2][2] = 1.0;
	//for ( int i = 0; i < 3; i++ )
	//{
	//	Fc_tau[i][i] += 1.0;
	//	Fc_t[i][i] += 1.0;
	//}
	
	// 2D plane-strain implementation detF
	detFc_t   =     Fc_t[0][0]*Fc_t[1][1] - Fc_t[0][1]*Fc_t[1][0];
	detFc_tau = Fc_tau[0][0]*Fc_tau[1][1] - Fc_tau[0][1]*Fc_tau[1][0];
	//
	// With the deformation gradient known at the element centroid
	// we are now able to implement the `F-bar' method later
	// ==================================================================
	
	int jj;
	double * svars  = SVARS; // current block of state variables
	// loop through integration points
	for ( int lquad = 0; lquad < nquad; ++lquad )
	{
		
		// ================================================
		// Obtain state variables from previous increment
        // ================================================
		if((*KINC <= 1)&&(*KSTEP == 1)) {
            
            // This is the first increment of the first step.
            //  Give initial conditions.
			onem(Fe_t);
	
        } else {
			
			// ===================================================
			// This is not the first increment; read old values.
            // get previous converged parameters
		    // ===================================================
			jj = 0;
			// elastic deformation gradient
            for ( int i = 0; i < 3; ++i ) {
				for ( int j = 0; j < 3; ++j ) {
					Fe_t[i][j] = svars[jj];
					jj++;
				}
			}
			
        }

		// ==================================================================
		// Obtain, and modify the deformation gradient at this integration
        //  point.  Modify the deformation gradienet for use in the `F-bar'
        //  method.  Also, take care of plane-strain.
		//
        onem(F_tau); // current
        onem(F_t);   // previously converged
        //zeroM2(F_tau);
		//zeroM2(F_t);
	
		for ( int i=0; i < ndm; ++i ) {
			for ( int j=0; j < ndm; ++j ) {
				for ( int k=0; k < nen; ++k ) {
					F_tau[i][j] += u[k][i]*dshp[lquad][k][j];
					F_t[i][j] += uOld[k][i]*dshp[lquad][k][j];
				}
			}
		}
		// modify for plane-strain
		F_tau[2][2] = 1.0;
		F_t[2][2] = 1.0;
		//for ( int i = 0; i < 3; i++ )
		//{
		//	F_tau[i][i] += 1.0;
		//	F_t[i][i] += 1.0;
		//}
		
        // Modify the deformation gradient for the `F-bar' method
        // only when using the 4 node fully integrated linear
        // element, do not use the `F-bar' method for any other element
		//
		// Note: LFLAGS[1]==1 means large-displacement analysis
        //
		if ((nen==4)&&(nquad==4)&&(LFLAGS[1]==1)) 
		{
			//  2D plane-strain implementation
            detF_t   =     F_t[0][0]*F_t[1][1] - F_t[0][1]*F_t[1][0];
			detF_tau = F_tau[0][0]*F_tau[1][1] - F_tau[0][1]*F_tau[1][0];
			double fac_t   = pow((detFc_t/detF_t),0.5);
			double fac_tau = pow((detFc_tau/detF_tau),0.5);
			
			// note below the factor does not multiply F33 (which remains at 1.0)
			for ( int i = 0; i < ndm; ++i ) {
				for ( int j = 0; j < ndm; ++j ) { 
					F_t[i][j]   = fac_t*F_t[i][j];
					F_tau[i][j] = fac_tau*F_tau[i][j];
					
				}
			}
			
        } else {
			cout << "error f-bar type\n";
			exit(1);
		}
		
		// Perform the constitutive time integration at this integ. point 
		// Input:
		//   PROPS  = material properties
		//   NPROPS = number of properties
		//   Fe_t   = previously converged elastic def gradient
		//   F_t    = previously converged total def gradient (F-bar previous)
		//   F_tau  = current total def gradient (F-bar current)
		//
		// Output:
		//   Fe_tau = current elastic def gradient
		//   T_tau    = current Cauchy stress
		//   SpTanMod = current tangent moduli
		//
		constitutive_update(PROPS,Fe_t,F_t,F_tau,
							Fe_tau,T_tau,SpTanMod);

		// Compute/update the displacement residual vector
        Smat[0] = T_tau[0][0];
		Smat[1] = T_tau[1][1];
		Smat[2] = T_tau[0][1];
		
		zeroM2(Bmat);
		for ( int a=0; a < nen; a++ )
		{
			Bmat[0][ndm*a  ] = dshpC[lquad][a][0]; // dNi/dx
			Bmat[1][ndm*a+1] = dshpC[lquad][a][1]; // dNi/dy
			Bmat[2][ndm*a  ] = dshpC[lquad][a][1]; // dNi/dy
			Bmat[2][ndm*a+1] = dshpC[lquad][a][0]; // dNi/dx
		}
		
		//BodyForceRes = zero
        // do kk=1,nNode
        //    BodyForceRes(1+nDim*(kk-1),1) = sh(kk)*body(1)
        //    BodyForceRes(2+nDim*(kk-1),1) = sh(kk)*body(2)
        // enddo
		
		// Swansea style
		for ( int a = 0; a < nen; a++ ) { // loop over nodes
			for ( int i = 0; i < ndm; ++i ) { // loop over degrees of freedom
				for ( int j = 0; j < 3; ++j ) { // loop over stress terms
					Ru[a][i] -= Bmat[j][ndm*a+i]*Smat[j]*detjC[lquad]*w[lquad];
				}
			} 
		} // end loop over nodes 		

		// Compute/update the displacement tangent matrix
        
        zeroM2(Gmat);
		for ( int a=0; a < nen; a++ )
		{
			Gmat[0][ndm*a  ] = dshpC[lquad][a][0]; // dNi/dx
			Gmat[1][ndm*a+1] = dshpC[lquad][a][0]; // dNi/dx
			Gmat[2][ndm*a  ] = dshpC[lquad][a][1]; // dNi/dy
			Gmat[3][ndm*a+1] = dshpC[lquad][a][1]; // dNi/dy
		}
		
		zeroM2(G0mat); // tangent moduli at centroid (note: only at centroid, index = 0)
		for ( int a=0; a < nen; a++ )
		{
			G0mat[0][ndm*a  ] = dshpC0[0][a][0]; // dNi/dx at centroid
			G0mat[1][ndm*a+1] = dshpC0[0][a][0]; // dNi/dx at centroid
			G0mat[2][ndm*a  ] = dshpC0[0][a][1]; // dNi/dy at centroid
			G0mat[3][ndm*a+1] = dshpC0[0][a][1]; // dNi/dy at centroid
		}
			
		zeroM2(Amat); // tangent moduli at current integration point
		
		Amat[0][0] = SpTanMod[0][0][0][0]; // 1111 -> 11 -> node1, dNi/dx (must jibe with Gmat above)
		Amat[0][1] = SpTanMod[0][0][1][0]; // 1121 -> 21 -> node2, dNi/dx
		Amat[0][2] = SpTanMod[0][0][0][1]; // 1112 -> 12 -> node1, dNi/dy
		Amat[0][3] = SpTanMod[0][0][1][1]; // 1122 -> 22 -> node2, dNi/dy
		
		Amat[1][0] = SpTanMod[1][0][0][0]; // 2111
		Amat[1][1] = SpTanMod[1][0][1][0]; // 2121
		Amat[1][2] = SpTanMod[1][0][0][1]; // 2112
		Amat[1][3] = SpTanMod[1][0][1][1]; // 2122
		
		Amat[2][0] = SpTanMod[0][1][0][0]; // 1211
		Amat[2][1] = SpTanMod[0][1][1][0]; // 1221
		Amat[2][2] = SpTanMod[0][1][0][1]; // 1212
		Amat[2][3] = SpTanMod[0][1][1][1]; // 1222
		
        Amat[3][0] = SpTanMod[1][1][0][0]; // 2211
		Amat[3][1] = SpTanMod[1][1][1][0]; // 2221
		Amat[3][2] = SpTanMod[1][1][0][1]; // 2212
		Amat[3][3] = SpTanMod[1][1][1][1]; // 2222
		
        zeroM2(Qmat); // tangent moduli at current integration point
		Qmat[0][0] = 0.5*( Amat[0][0]+Amat[0][3] ) - 0.5*T_tau[0][0];
		Qmat[1][0] = 0.5*( Amat[1][0]+Amat[1][3] ) - 0.5*T_tau[0][1];
		Qmat[2][0] = 0.5*( Amat[2][0]+Amat[2][3] ) - 0.5*T_tau[0][1];
		Qmat[3][0] = 0.5*( Amat[3][0]+Amat[3][3] ) - 0.5*T_tau[1][1];
		
		Qmat[0][3] = 0.5*( Amat[0][0]+Amat[0][3] ) - 0.5*T_tau[0][0];
		Qmat[1][3] = 0.5*( Amat[1][0]+Amat[1][3] ) - 0.5*T_tau[0][1];
		Qmat[2][3] = 0.5*( Amat[2][0]+Amat[2][3] ) - 0.5*T_tau[0][1];
		Qmat[3][3] = 0.5*( Amat[3][0]+Amat[3][3] ) - 0.5*T_tau[1][1];

		
		
		
	
		// Swansea style
		assert(nsdm == 4);
		/*for (int a=0; a<nen; a++) { // loop over nodes
			for(int i=0; i<ndf; i++){ // loop over dofs of node
				for(int b=0; b<nen; b++) { // loop over nodes
					for(int j=0; j<ndf; j++) { // loop over dofs of node					
						for(int I=0; I<nsdm; I++) { // loop over nsdm == ndm*ndf (tangent)
							for(int J=0; J<nsdm; J++) { // loop over nsdm == ndm*ndf (tangent)
								Kuu[a*ndf+i][b*ndf+j]+= Gmat[I][ndf*a+i]*Amat[I][J]*Gmat[J][ndf*b+j]*detjC[lquad]*w[lquad];														
							}
						}
					}
				}
			}
		}
		
		for (int a=0; a<nen; a++) { // loop over nodes
			for(int i=0; i<ndf; i++){ // loop over dofs of node
				for(int b=0; b<nen; b++) { // loop over nodes
					for(int j=0; j<ndf; j++) { // loop over dofs of node					
						for(int I=0; I<nsdm; I++) { // loop over nsdm == ndm*ndf (tangent)
							for(int J=0; J<nsdm; J++) { // loop over nsdm == ndm*ndf (tangent)
								Kuu[a*ndf+i][b*ndf+j]+= Gmat[I][ndf*a+i]*Qmat[I][J]*( G0mat[J][ndf*b+j]-Gmat[J][ndf*b+j] )*detjC[lquad]*w[lquad];
							}
						}
					}
				}
			}
		}*/
		
		for(int i=0; i <nen*ndf; i++){ // loop over dofs of node
			for(int j=0; j <nen*ndf; j++) { // loop over dofs of node					
				for(int II=0; II <nsdm; II++) { // loop over nsdm == ndm*ndf (tangent)
					for(int JJ=0; JJ <nsdm; JJ++) { // loop over nsdm == ndm*ndf (tangent)
						Kuu[i][j] += detjC[lquad]*w[lquad]*Gmat[II][i]*(
							Amat[II][JJ]*Gmat[JJ][j]
						+   Qmat[II][JJ]*(G0mat[JJ][j]-Gmat[JJ][j])
							);
							
					}
				}
			}
		}
		
		

        // ================================================
		// Store current state (Abaqus will 
		// automatically perform update at convergence)
		// ================================================
		jj = 0;
		
		// elastic deformation gradient
        for ( int i = 0; i < 3; ++i ) {
			for ( int j = 0; j < 3; ++j ) {
				svars[jj] = Fe_tau[i][j];
				jj++;
			}
		}
		
		// ================================================
		// Take strides
        // ================================================
		svars += nlSdv; // setup for the next intPt

		
	} // end loop over integration points

	/*cout << "=====================================================" << endl;
    cout << "Kuu at step = " << *KSTEP << " inc = " << *KINC << endl;
	for ( int i = 0; i < 8; ++i )
	{
		for ( int j = 0; j < 8; ++j )
		{
			cout.width(20); 
			cout.precision(18);
			cout << std::right << Kuu[i][j] << " ";
		
		}
		cout << endl;
    }*/
	
	// Assemble the element level residual
	for ( int a = 0; a < nen; a++ ) {
		for(int i = 0; i < ndf; i++ ) { // loop over dofs of node
			RHS[a*ndf+i] = Ru[a][i];
		}
	}
	
	// Assemble the element level tangent matrix
	// note: fortran is column major i.e.
	// amatrx = 
	// 0  4   8  12 
	// 1  5   9  13
	// 2  6  10  14
	// 3  7  11  15
	// Kuu = 
	// (00) (01) (02) (03)
	// (10) (11) (12) (13)
	// (20) (21) (22) (23)
	// (30) (31) (32) (33)
	// Therefore:
	// amatrx[row + col*(nen*ndf)] = Kuu[row][col]
	int row,col;
	for (int a=0; a<nen; a++) { // loop over nodes
		for(int i=0; i<ndf; i++){ // loop over dofs of node
			for(int b=0; b<nen; b++) { // loop over nodes
				for(int j=0; j<ndf; j++) { // loop over dofs of node					
					row = a*ndf+i;
					col = b*ndf+j;
					FORT_ACCESS(AMATRX,row,col,nen*ndf)=Kuu[row][col];							
				}
			}
		}
	}
				
	/*for ( int i = 0; i < nen*ndf; ++i ) {
		for ( int j = 0; j < nen*ndf; ++j ) {
			cout.width(10); 
			cout << std::right << Kuu[i][j] << " ";
		}
		cout << endl;
	}
	cin.get();*/
}

void constitutive_update(double * PROPS, Tensor2 & Fe_t, Tensor2 & F_t, Tensor2 & F_tau,
						 Tensor2 & Fe_tau, Tensor2 & T_tau, Tensor4 & SpTanMod)
{
	Tensor2  Iden(boost::extents[3][3]);       // identity
	
	Tensor2  F_t_inv(boost::extents[3][3]);    // inverse of previous converged def gradient
	Tensor2  finc(boost::extents[3][3]);       // incremental def gradient
	
	Tensor2  FeTr(boost::extents[3][3]);      // TRIAL elastic def gradient
	Tensor2  VeTr(boost::extents[3][3]);      // TRIAL elastic left strength
	VectorD  VePrTr(boost::extents[3]);       // eigenvalues of VeTr
	Tensor2  VeTr_inv(boost::extents[3][3]);  // inverse of TRIAL elastic left strength
	Tensor2  ReTr(boost::extents[3][3]);      // TRIAL elastic rotation
	
	Tensor2  Ve_tau(boost::extents[3][3]);      // stretch part of current ELASTIC def gradient
	VectorD  VePr_tau(boost::extents[3]);       // eigenvalues of Ve_tau
	Tensor2  Re_tau(boost::extents[3][3]);      // rotation part of Ve_tau
	
	Tensor2  epsE_tau(boost::extents[3][3]);      // elastic log strain at time tau
	VectorD  epsEPr_tau(boost::extents[3]);       // eigenvalues of epsE_tau
	MatrixD2 Nvec_epsE_tau(boost::extents[3][3]); // eigenvectors of epsE_tau
	
	//Tensor2  Be_t(boost::extents[3][3]);       // left Cauchy-Green tensor at previous converged time, Be_t
	//VectorD  BePr_t(boost::extents[3]);        // eigenvalues of Be_t
	
	Tensor2  BeTr(boost::extents[3][3]);        // TRIAL left Cauchy-Green tensor, BeTr
	VectorD  BeTrPr(boost::extents[3]);         // eigenvalues of BeTr
	MatrixD2 Nvec_BeTr(boost::extents[3][3]);   // eigenvectors of BeTr
	
	Tensor2  epsEtr(boost::extents[3][3]);      // TRIAL elastic log strain
	VectorD  epsEtrPr(boost::extents[3]);       // eigenvalues of TRIAL elastic log strain
	
	MatrixD2 Dalg(boost::extents[9][9]);         // small-strain CTO (size 9 x 9)
	Tensor2  kirSig(boost::extents[3][3]);       // Kirchhoff stress tensor (output from small-strain update)
	
	Tensor4  LTensor(boost::extents[3][3][3][3]); // LTensor 
	Tensor4  BTensor(boost::extents[3][3][3][3]); // BTensor 
	Tensor4  STensor(boost::extents[3][3][3][3]); // LTensor 
	
	MatrixD2 LMat(boost::extents[9][9]);         // LMat converted from LTensor
	MatrixD2 BMat(boost::extents[9][9]);         // BMat converted from BTensor 
	MatrixD2 SMat(boost::extents[9][9]);         // SMat converted from STensor 
	
	MatrixD2 DLB(boost::extents[9][9]);          // 1/(2*det(F))*Dalg*L*B 
	
	MatrixD2 AMat(boost::extents[9][9]);         // spatial tangent in 9x9 matrix form 
	
	double detF_tau;
	
	// identity
	onem(Iden);
	
	// inverse of previous converged def gradient
	matlib_inverse(F_t,F_t_inv);
	
	// compute incremental deformation gradient
	matlib_matmult(F_tau,F_t_inv,finc);
	
	// ===========================================================================
	// compute FeTr = finc*Fe_t
	// ===========================================================================
	matlib_matmult(finc,Fe_t,FeTr);
	
	// ===========================================================================
	// compute trial left Cauchy-Green tensor BeTr = FeTr * FeTr^T
	// ===========================================================================
	zeroT2(BeTr);
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			for ( int k = 0;  k < 3; ++k ) {
				BeTr[i][j] += FeTr[i][k]*FeTr[j][k];
			}
		}
	}
	
	// ===========================================================================
	// compute trial elastic (half) log strain tensor
	// ===========================================================================
	// get eigenvalues and eigenvectors of Be_TR
	solve_eigenvalue_3x3(BeTr, BeTrPr, Nvec_BeTr);

	for ( int i = 0; i < 3; ++i ) {
		epsEtrPr[i] = 0.5*log(BeTrPr[i]); // compute principal (half) log strains
		VePrTr[i] = sqrt(BeTrPr[i]);      // eigenvalues of VeTr
	}
	
    // compute epsEtr
	zeroT2(epsEtr);
	for ( int k = 0; k < 3; ++k ) {
		for ( int i = 0; i < 3; ++i ) {
			for ( int j = 0; j < 3; ++j ) {
				epsEtr[i][j] += epsEtrPr[k]*Nvec_BeTr[k][i]*Nvec_BeTr[k][j];
			}		
		}
	}
	
	// ===========================================================================
	// compute VeTr and ReTr parts of FeTr
	// ===========================================================================
	// compute TRIAL left stretch VeTr
	zeroT2(VeTr);
	for ( int k = 0; k < 3; ++k ) {
		for ( int i = 0; i < 3; ++i ) {
			for ( int j = 0; j < 3; ++j ) {
				VeTr[i][j] += VePrTr[k]*Nvec_BeTr[k][i]*Nvec_BeTr[k][j];
			}		
		}
	}
	
	// inverse of VeTr
	matlib_inverse(VeTr,VeTr_inv);
	
	// TRIAL elastic rotation
	matlib_matmult(VeTr_inv,FeTr,ReTr);
	
	// =======================================================================
	// small-strain return map to get small-strain consistent
    // tangent (Dalg), Kirchhoff stress, and elastic small-strain tensor
    // =======================================================================
	VMconst(PROPS,epsEtr,
			Dalg,kirSig,epsE_tau);
    
	// =======================================================================
	// convert epsE_tau to Fe_tau
	// =======================================================================
	// Re_tau = ReTR
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			Re_tau[i][j] = ReTr[i][j];
		}		
	}
	
	// Ve_tau = exp(epsE_tau)
	// get eigenvalues and eigenvectors of epsE_tau
	solve_eigenvalue_3x3(epsE_tau, epsEPr_tau, Nvec_epsE_tau);
	
	// eigenvalues of epsE_tau
	for ( int i = 0; i < 3; ++i ) {
		VePr_tau[i] = exp(epsEPr_tau[i]);      // eigenvalues of Ve_tau
	}
	
	zeroT2(Ve_tau);
	for ( int k = 0; k < 3; ++k ) {
		for ( int i = 0; i < 3; ++i ) {
			for ( int j = 0; j < 3; ++j ) {
				Ve_tau[i][j] += VePr_tau[k]*Nvec_epsE_tau[k][i]*Nvec_epsE_tau[k][j];
			}		
		}
	}
	
	
	// compute Fe_tau = Ve_tau * Re_tau
	matlib_matmult(Ve_tau,Re_tau,Fe_tau);
	
	/*// check Re_tau
	Tensor2 ReChk(boost::extents[3][3]);
	Tensor2 Ve_tau_inv(boost::extents[3][3]);  // inverse of Ve_tau
	matlib_inverse(Ve_tau,Ve_tau_inv);
	matlib_matmult(Ve_tau_inv,Fe_tau,ReChk);
	cout << "===============================\n";
	cout << "Re_tau = " << endl;
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			cout << Re_tau[i][j] << " ";
		}
		cout << endl;
	}
	
	cout << "ReChk = " << endl;
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			cout << ReChk[i][j] << " ";
		}
		cout << endl;
	}*/
	
	// ===========================================================================
	// compute J = determinant of F_tau
	// note: F_tau is the trial deformation gradient i.e. F_tau = finc*F_t
    detF_tau = matlib_determinant(F_tau);
	
	// scale Dalg by 1/(2*det(F_tau))
	for ( int i = 0; i < 9; ++i ) {
		for ( int j = 0; j < 9; ++j ) {
			Dalg[i][j] = Dalg[i][j]/(2.0*detF_tau);
		}
	}
		
	// convert Kirchhoff stress to Cauchy stress
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			T_tau[i][j] = kirSig[i][j]/detF_tau;
		}
	}
	
	// compute LMat
	parDerGen(BeTr,BeTrPr,Nvec_BeTr,LTensor);
	convertTensor4ToMat9(LTensor,LMat);
	
	// compute B_ijkl
	// B_ijkl = delta_ik B^e,tr_jl + delta_jk B^e,tr_il (see Eq. 14.102, pg. 598 in complas peric)
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			for ( int k = 0; k < 3; ++k ) {
				for ( int l = 0; l < 3; ++l ) {
					BTensor[i][j][k][l] = Iden[i][k]*BeTr[j][l] + Iden[j][k]*BeTr[i][l];
				}
			}
		}
	}
	convertTensor4ToMat9(BTensor, BMat);
	
	// compute S_ijkl = sigma_il * delta_jk
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			for ( int k = 0; k < 3; ++k ) {
				for ( int l = 0; l < 3; ++l ) {
					STensor[i][j][k][l] = T_tau[i][l]*Iden[j][k];
				}
			}
		}
	}
	convertTensor4ToMat9(STensor, SMat);
	
	// compute 1/(2*det(F))*Dalg*L*B
	zeroM2(DLB);
	for ( int i = 0; i < 9; ++i ) {
		for ( int j = 0; j < 9; ++j ) {
			for ( int m = 0; m < 9; m++ ) {
				for ( int n = 0; n < 9; n++ ) {
					DLB[i][j] += Dalg[i][m]*LMat[m][n]*BMat[n][j];
				}
			}
		}
	}
	
	// spatial tangent in 9x9 matrix form
	//a=1/(2*det(F))*Dalg*L*B-S;
	for ( int i = 0; i < 9; ++i ) {
		for ( int j = 0; j < 9; ++j ) {
			AMat[i][j] = DLB[i][j] - SMat[i][j];
		}
	}
	
	// convert spatial tangent to tensor form
	convertMat9ToTensor4(AMat,SpTanMod);
	
}

void parDerGen(Tensor2 & X, VectorD & eP, MatrixD2 & eV, Tensor4 & L)
{
	double tol=1.e-12; 
	VectorD yP(boost::extents[3]);       
	VectorD ydash(boost::extents[3]);       
	Tensor2 Iden(boost::extents[3][3]);
	Tensor4 Is(boost::extents[3][3][3][3]);
	
	onem(Iden);
	
    for (int i=0; i < 3; ++i ) {
		for (int j=0; j < 3; ++j ) {
			for (int k=0; k < 3; ++k ) {
				for (int l=0; l < 3; ++l ) {
					Is[i][j][k][l] = 0.5*(Iden[i][k]*Iden[j][l]+Iden[i][l]*Iden[j][k]);
				}
			}
		}
	}
	
	for ( int i = 0; i < 3; ++i ) {
		yP[i]    = log(eP[i]);
		ydash[i] = 1.0/eP[i];
	}
    
    
    if (abs(eP[0])<tol && abs(eP[1])<tol && abs(eP[2])<tol) { // all 3 eigenvalues are ZERO
    
		//L=Is;
		for (int i=0; i < 3; ++i ) {
			for (int j=0; j < 3; ++j ) {
				for (int k=0; k < 3; ++k ) {
					for (int l=0; l < 3; ++l ) {
						L[i][j][k][l]=Is[i][j][k][l];
					}
				}
			}
		}
		
	} else if ( abs(eP[0]-eP[1])<tol && abs(eP[0]-eP[2])<tol ) { // all 3 eigenvalues are equal but NOT ZERO
    
		// all eigenvalues are (approx.) the SAME
        // d_ab = d(ln \lam_a)/d(\lam_b) = 1/lam_a * delta_ab
        // and alternative form of outlined by Ogden,
        // 0.5*( ln \lam_a - ln \lam_b )/( \lam_a - \lam_b ) 
        // is replaced by 0.5*d( ln \lam_a - ln \lam_b )/d( \lam_a  ) 
        // i.e. from l'Hospital rule
        // = 0.5*( 1 - delta_ab )/\lam_a = 0.5/\lam_a for a != b	
        // take derivative directly on G matrix
        // note: ydash[0]=1/BePr[0];    
		//L=ydash[0]*Is;
		for (int i=0; i < 3; ++i ) {
			for (int j=0; j < 3; ++j ) {
				for (int k=0; k < 3; ++k ) {
					for (int l=0; l < 3; ++l ) {
						L[i][j][k][l]=ydash[0]*Is[i][j][k][l];
					}
				}
			}
		}
	
		
	} else if ( abs(eP[0]-eP[1])<tol || abs(eP[1]-eP[2])<tol || abs(eP[0]-eP[2])<tol ) { // 2 eigenvalues are equal
        
		double xa,xc;
		double ya,yc;
		double yda,ydc;
        if ( abs(eP[0]-eP[1])<tol ) {
            // a=2 != b=0,c=1
            xa=eP[2]; 
            xc=eP[0]; 
            ya=yP[2]; 
            yc=yP[0]; 
            yda=ydash[2]; 
            ydc=ydash[0]; 
		} else if ( abs(eP[1]-eP[2])<tol ) {
            xa=eP[0]; 
            xc=eP[1]; 
            ya=yP[0]; 
            yc=yP[1]; 
            yda=ydash[0]; 
            ydc=ydash[1]; 
		} else {                       
            xa=eP[1]; 
            xc=eP[0]; 
            ya=yP[1]; 
            yc=yP[0]; 
            yda=ydash[1]; 
            ydc=ydash[0]; 
		}
		//
		// ya = log(lambda_a)
		// xa = lambda_a
		// dy_a/dx_b = 0 for a not eq. to b
		// dy_a/dx_a = 1/lambda_a for a == b
		// apply this to Eq. A.47 in complas peric
		//
		// yda = dy_a/dx_a
		// ydb = dy_b/dx_b
		// ydc = dy_c/dx_c
		// and cross terms (e.g. dy_a/dx_b) are zero
		//
		Tensor4 dX2dX(boost::extents[3][3][3][3]);
		Tensor4 XxX(boost::extents[3][3][3][3]);
		Tensor4 XxI(boost::extents[3][3][3][3]);
		Tensor4 IxX(boost::extents[3][3][3][3]);
		Tensor4 IxI(boost::extents[3][3][3][3]);
		VectorD s(boost::extents[6]);
		
		for (int i=0; i < 3; ++i ) {
			for (int j=0; j < 3; ++j ) {
				for (int k=0; k < 3; ++k ) {
					for (int l=0; l < 3; ++l ) {
						dX2dX[i][j][k][l] = 0.5*(Iden[i][k]*X[l][j]+Iden[i][l]*X[k][j]+Iden[j][l]*X[i][k]+Iden[k][j]*X[i][l]);
						XxX[i][j][k][l] = X[i][j]*X[k][l];
						XxI[i][j][k][l] = X[i][j]*Iden[k][l];
						IxX[i][j][k][l] = Iden[i][j]*X[k][l];
						IxI[i][j][k][l] = Iden[i][j]*Iden[k][l];
					}
				}
			}
		}
		
        

        // see Eq. A.53
		double xaxc2 = (xa-xc)*(xa-xc);
		double xaxc3 = xaxc2*(xa-xc);
		
        s[0]=(ya-yc)/xaxc2-ydc/(xa-xc);  
        s[1]=2.0*xc*(ya-yc)/xaxc2-(xa+xc)/(xa-xc)*ydc;
        s[2]=2.0*(ya-yc)/xaxc3-(yda+ydc)/xaxc2; 
        s[3]=xc*s[2];    // eq. to s4 and s5 in A.47
        s[4]=s[3];       // eq. to s4 and s5 in A.47
        s[5]=xc*xc*s[2]; // eq. to s6 in A.47

        //L = s[0]*dX2dX - s[1]*Is - s[2]*XxX + s[3]*XxI + s[4]*IxX - s[5]*IxI;
		for (int i=0; i < 3; ++i ) {
			for (int j=0; j < 3; ++j ) {
				for (int k=0; k < 3; ++k ) {
					for (int l=0; l < 3; ++l ) {
						L[i][j][k][l] = s[0]*dX2dX[i][j][k][l]
									    - s[1]*Is[i][j][k][l] 
									    - s[2]*XxX[i][j][k][l]
										+ s[3]*XxI[i][j][k][l]
										+ s[4]*IxX[i][j][k][l]
										- s[5]*IxI[i][j][k][l];
		
					}
				}
			}
		}
	
		
	} else { //distinct eigenvalues
        
		Tensor4 dX2dX(boost::extents[3][3][3][3]);
		
		for (int i=0; i < 3; ++i ) {
			for (int j=0; j < 3; ++j ) {
				for (int k=0; k < 3; ++k ) {
					for (int l=0; l < 3; ++l ) {
						dX2dX[i][j][k][l] = 0.5*(Iden[i][k]*X[l][j]+Iden[i][l]*X[k][j]+Iden[j][l]*X[i][k]+Iden[k][j]*X[i][l]);
					}
				}
			}
		}
		
		
        // eigenprojections
        
        MatrixD3 E(boost::extents[3][3][3]);
		zeroM3(E);
		for ( int k = 0; k < 3; ++k ) {
			for ( int i = 0; i < 3; ++i ) {
				for ( int j = 0; j < 3; ++j ) {
					E[k][i][j] += eV[k][i]*eV[k][j];
				}		
			}
		}
	
		Tensor4 EaxEa(boost::extents[3][3][3][3]);
		Tensor4 EbxEb(boost::extents[3][3][3][3]);
		Tensor4 EcxEc(boost::extents[3][3][3][3]);
		
        MatrixD2 permut(boost::extents[3][3]);
		
		permut[0][0] = 0;
		permut[0][1] = 1;
		permut[0][2] = 2;
		
		permut[1][0] = 1;
		permut[1][1] = 2;
		permut[1][2] = 0;
		
		permut[2][0] = 2;
		permut[2][1] = 0;
		permut[2][2] = 1;
		
        
		int i1,i2,i3;
		double xa,xb,xc;
		double ya;
		double yda;
        
		zeroT4(L);
		for ( int p = 0; p <3; ++p ) 
		{
            
            i1 = permut[p][0];
            i2 = permut[p][1];
            i3 = permut[p][2];
			
			for (int i=0; i < 3; ++i ) {
				for (int j=0; j < 3; ++j ) {
					for (int k=0; k < 3; ++k ) {
						for (int l=0; l < 3; ++l ) {
							EaxEa[i][j][k][l] = E[i1][i][j]*E[i1][k][l];
							EbxEb[i][j][k][l] = E[i2][i][j]*E[i2][k][l];
							EcxEc[i][j][k][l] = E[i3][i][j]*E[i3][k][l];
						}
					}
				}
			}
					
            xa=eP[i1]; 
            xb=eP[i2]; 
            xc=eP[i3]; 
            ya=yP[i1]; 
            yda=ydash[i1]; 

			// L += ya*( dX2dX
			// 		  - (xb+xc)*Is
			// 		  - ((xa-xb)+(xa-xc))*EaxEa
			// 		  - (xb-xc)*(EbxEb - EcxEc) 
			// 		   )/((xa-xb)*(xa-xc)) + yda*EaxEa;

			for (int i=0; i < 3; ++i ) {
				for (int j=0; j < 3; ++j ) {
					for (int k=0; k < 3; ++k ) {
						for (int l=0; l < 3; ++l ) {
							L[i][j][k][l] += ya*( dX2dX[i][j][k][l] 
												- (xb+xc)*Is[i][j][k][l]
												- ((xa-xb)+(xa-xc))*EaxEa[i][j][k][l]
												- (xb-xc)*(EbxEb[i][j][k][l] - EcxEc[i][j][k][l]) 
												)/((xa-xb)*(xa-xc)) + yda*EaxEa[i][j][k][l];

						}
					}
				}
			}
			

		}
		
	} // end cases
	
	
}

void VMconst(double * PROPS, Tensor2 & epsEtr,
			 MatrixD2 & Dalg, Tensor2 & sig_tau, Tensor2 & epsE_tau)
{
	
	double E,nu,sigmay;
	
	// Obtain material properties
    E      = PROPS[0];
    nu     = PROPS[1];
	sigmay = PROPS[2];
	
	// ===============================================================================================	
	double tol = 1.e-12;
	double Kmod,Gmod;
	double SdevNorm,trSig,j2,f;
	
	EgVecD bm1; bm1.resize(9); bm1.setZero();
	
	EgVecD epsE9(9); 
	EgVecD epsEtr9(9); 
	EgVecD Sig9(9);
	EgVecD Sdev9(9); 
	EgVecD Nvec(9); 
	
	EgMatD2 Dalg9(9,9); Dalg9.setZero();
	EgMatD2 De9(9,9);   De9.setZero();
	EgMatD2 Ce9(9,9);   Ce9.setZero();
	EgMatD2 Id9(9,9);   Id9.setZero();
	EgMatD2 Idev9(9,9);  Idev9.setZero();
	EgMatD2 Pmat(9,6); Pmat.setZero(); // matrix to convert between 9- and 6-element stuff
	
	// identity tensor size-9 vector form
	bm1(0) = 1.0;
	bm1(1) = 1.0;
	bm1(2) = 1.0;
	
	Id9  = EgMatD2::Identity(9,9);
	Idev9 = Id9 - bm1*bm1.transpose()/3.0;
	
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
	
	// Trial Cauchy stress
	Sig9 = De9*epsEtr9;
	
	// Trial deviatoric stress vector
	trSig = (Sig9(0)+Sig9(1)+Sig9(2))/3.0;
	Sdev9 = Sig9-trSig*bm1;
	
	// Compute J2=1/2*s_ij*s_ij
	SdevNorm = Sdev9.norm();
	j2 = 0.5*SdevNorm*SdevNorm;
	
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
		maxit = 3;
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
		//
		// |s| = sqrt(2)*sqrt(j2)
		//
		//df=sqrt(3)/(2*fc*sqrt(j2))*dj2;
		//  =sqrt(3)/(2*fc*|s|/sqrt(2))*dj2
		//  =sqrt(3)/(sqrt(2)*fc*|s|)*dj2
		//  =sqrt(3)/(sqrt(2)*fc)*dj2/|s|
		Nvec = Sdev9/SdevNorm;
		df  = sqrt(3.0/2.0)/sigmay*Nvec;
		ddf = sqrt(3)/(2.0*sigmay*sqrt(j2))*( Idev9 - Nvec*Nvec.transpose() ); 
   
       
		// newton iteration
		for ( int itnum = 1; itnum <= maxit; ++itnum )
		{
			//cout << "\tnewton iter = " << itnum << " -------------------- " << endl;
			
			// form Jacobian
			A.setZero();
			A.topLeftCorner(9,9) = Id9+dgam*ddf*De9;
			A.block(9,0,1,9) = df.transpose()*De9;
			A.block(0,9,9,1) = df;
				
			Ainv = A.inverse();
			
			dx = -Ainv * resid;
			
			for ( int i = 0; i <=8; ++i ) {
				epsE9(i) += dx(i); // elems 1 through 9
			}			
			dgam  += dx(9);      // elem 10
			Sig9 = De9*epsE9;
			
			// deviatoric stress vector
			trSig = (Sig9(0)+Sig9(1)+Sig9(2))/3.0;
			Sdev9 = Sig9-trSig*bm1;
	
			// Compute J2=1/2*s_ij*s_ij
			SdevNorm = Sdev9.norm();
			j2 = 0.5*SdevNorm*SdevNorm;
	
			// recompute terms in jacobian matrix for next iteration
			Nvec = Sdev9/SdevNorm;
			df  = sqrt(3.0/2.0)/sigmay*Nvec;
			ddf = sqrt(3.0)/(2.0*sigmay*sqrt(j2))*( Idev9 - Nvec*Nvec.transpose() ); 
   
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
			
			// consistent tangent operator
			A.setZero();
			A.topLeftCorner(9,9) = Ce9 + dgam*ddf;
			A.block(9,0,1,9) = df.transpose();
			A.block(0,9,9,1) = df;
							
			EgMatD2 B = A.inverse();
			Dalg9 = B.topLeftCorner(9,9); 
			
			//cout << "Dalg6 = " << endl << Dalg6 << endl;
			
		} else {
			cout << "no convergence in return map after " << maxit << " iterations " << endl;
			exit(1);	
		}
		
	} else {
		
		// ===============================
		// ELASTIC UPDATE
		// ===============================
		cout << "Elastic update\n";
		Dalg9 = De9;
		

	} // end check yield	
	
	//cout << "*****************************************************\n";
	//cout << "epsE9 = " << epsE9 << endl;
	//cout << "*****************************************************\n";
	
	// convert strains from 9 to boost:multiarray tensor form
	epsE_tau[0][0] = epsE9(0); // xx
	epsE_tau[1][1] = epsE9(1); // yy
	epsE_tau[2][2] = epsE9(2); // zz
	epsE_tau[0][1] = epsE9(3); // xy
	epsE_tau[1][0] = epsE9(4); // yx
	epsE_tau[1][2] = epsE9(5); // yz
	epsE_tau[2][1] = epsE9(6); // zy
	epsE_tau[2][0] = epsE9(7); // zx
	epsE_tau[0][2] = epsE9(8); // xz
	
	// convert stresses from 9 to output form
	sig_tau[0][0] = Sig9(0); // xx
	sig_tau[1][1] = Sig9(1); // yy
	sig_tau[2][2] = Sig9(2); // zz
	sig_tau[0][1] = Sig9(3); // xy
	sig_tau[1][0] = Sig9(4); // yx
	sig_tau[1][2] = Sig9(5); // yz
	sig_tau[2][1] = Sig9(6); // zy
	sig_tau[2][0] = Sig9(7); // zx
	sig_tau[0][2] = Sig9(8); // xz
	
	// assign Eigen++ Dalg9 to Boost Dalg
	for ( int i = 0; i < 9; ++i ) {
		for ( int j = 0; j < 9; ++j ) {
			Dalg[i][j] = Dalg9(i,j);
		}
	}
}

// -------------------------------------------------
// CONVENIENCE FUNCTIONS
// -------------------------------------------------
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
double matlib_determinant(Tensor2 & A)
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
double matlib_inverse(Tensor2 & A, Tensor2 & Ainv)
{
	double det;

	det = matlib_determinant(A);
	if (fabs(det) < 1.e-20) return 0.e0;

	
	// 0 1 2 = [0][0] [0][1] [0][2]
	// 3 4 5 = [1][0] [1][1] [1][2]
	// 6 7 8 = [2][0] [2][1] [2][2]

	Ainv[0][0] = ( A[1][1]*A[2][2]-A[1][2]*A[2][1])/det;
	Ainv[0][1] = (-A[0][1]*A[2][2]+A[0][2]*A[2][1])/det;
	Ainv[0][2] = ( A[0][1]*A[1][2]-A[0][2]*A[1][1])/det;
	Ainv[1][0] = (-A[1][0]*A[2][2]+A[1][2]*A[2][0])/det;
	Ainv[1][1] = ( A[0][0]*A[2][2]-A[0][2]*A[2][0])/det;
	Ainv[1][2] = (-A[0][0]*A[1][2]+A[0][2]*A[1][0])/det;
	Ainv[2][0] = ( A[1][0]*A[2][1]-A[1][1]*A[2][0])/det;
	Ainv[2][1] = (-A[0][0]*A[2][1]+A[0][1]*A[2][0])/det;
	Ainv[2][2] = ( A[0][0]*A[1][1]-A[0][1]*A[1][0])/det;

	return det;
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
	}
	
}
void convertTensor4ToMat9(Tensor4 & T, MatrixD2 & A)
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

void convertMat9ToTensor4(MatrixD2 & A, Tensor4 & T)
{
	T[0][0][0][0] = A[0][0];
	T[0][0][1][1] = A[0][1];
	T[0][0][2][2] = A[0][2];
	T[0][0][0][1] = A[0][3];
	T[0][0][1][0] = A[0][4];
	T[0][0][1][2] = A[0][5];
	T[0][0][2][1] = A[0][6];
	T[0][0][2][0] = A[0][7];
	T[0][0][0][2] = A[0][8];
	             
	T[1][1][0][0] = A[1][0];
	T[1][1][1][1] = A[1][1];
	T[1][1][2][2] = A[1][2];
	T[1][1][0][1] = A[1][3];
	T[1][1][1][0] = A[1][4];
	T[1][1][1][2] = A[1][5];
	T[1][1][2][1] = A[1][6];
	T[1][1][2][0] = A[1][7];
	T[1][1][0][2] = A[1][8];
	             
	T[2][2][0][0] = A[2][0];
	T[2][2][1][1] = A[2][1];
	T[2][2][2][2] = A[2][2];
	T[2][2][0][1] = A[2][3];
	T[2][2][1][0] = A[2][4];
	T[2][2][1][2] = A[2][5];
	T[2][2][2][1] = A[2][6];
	T[2][2][2][0] = A[2][7];
	T[2][2][0][2] = A[2][8];
	             
	T[0][1][0][0] = A[3][0];
	T[0][1][1][1] = A[3][1];
	T[0][1][2][2] = A[3][2];
	T[0][1][0][1] = A[3][3];
	T[0][1][1][0] = A[3][4];
	T[0][1][1][2] = A[3][5];
	T[0][1][2][1] = A[3][6];
	T[0][1][2][0] = A[3][7];
	T[0][1][0][2] = A[3][8];
	             
	T[1][0][0][0] = A[4][0];
	T[1][0][1][1] = A[4][1];
	T[1][0][2][2] = A[4][2];
	T[1][0][0][1] = A[4][3];
	T[1][0][1][0] = A[4][4];
	T[1][0][1][2] = A[4][5];
	T[1][0][2][1] = A[4][6];
	T[1][0][2][0] = A[4][7];
	T[1][0][0][2] = A[4][8];
	             
	T[1][2][0][0] = A[5][0];
	T[1][2][1][1] = A[5][1];
	T[1][2][2][2] = A[5][2];
	T[1][2][0][1] = A[5][3];
	T[1][2][1][0] = A[5][4];
	T[1][2][1][2] = A[5][5];
	T[1][2][2][1] = A[5][6];
	T[1][2][2][0] = A[5][7];
	T[1][2][0][2] = A[5][8];
	             
	T[2][1][0][0] = A[6][0];
	T[2][1][1][1] = A[6][1];
	T[2][1][2][2] = A[6][2];
	T[2][1][0][1] = A[6][3];
	T[2][1][1][0] = A[6][4];
	T[2][1][1][2] = A[6][5];
	T[2][1][2][1] = A[6][6];
	T[2][1][2][0] = A[6][7];
	T[2][1][0][2] = A[6][8];
	             
	T[2][0][0][0] = A[7][0];
	T[2][0][1][1] = A[7][1];
	T[2][0][2][2] = A[7][2];
	T[2][0][0][1] = A[7][3];
	T[2][0][1][0] = A[7][4];
	T[2][0][1][2] = A[7][5];
	T[2][0][2][1] = A[7][6];
	T[2][0][2][0] = A[7][7];
	T[2][0][0][2] = A[7][8];
	             
	T[0][2][0][0] = A[8][0];
	T[0][2][1][1] = A[8][1];
	T[0][2][2][2] = A[8][2];
	T[0][2][0][1] = A[8][3];
	T[0][2][1][0] = A[8][4];
	T[0][2][1][2] = A[8][5];
	T[0][2][2][1] = A[8][6];
	T[0][2][2][0] = A[8][7];
	T[0][2][0][2] = A[8][8];
	
	
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

