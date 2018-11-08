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
!  *User Element,Nodes=4,Type=U1,Iproperties=?,Properties=2,Coordinates=2,Variables=9,Unsymm
!  1,2
!
! Note: 9 state variables are used in this element to store previously 
! converged elastic deformation gradient (Fe-bar)
!
! No integer properties (jprops) used for this element - we only use nInt = 4 (full integration) 
!
! Material Properties Vector
! --------------------------------------------------------------
! lambda  = props(1)  ! K
! mu      = props(2)  ! G
!

The C++ umat must be compiled externally using a compiler supported by ABAQUS, producing an object file .obj. The compiler command is

cl /EHsc /c /Fo uel4nodeElastic.cpp /I"C:\cygwin\home\KWL\libraries\eigen-3.2.1" /I"C:\Users\KWL\Desktop\Research\libraries\boost_1_56_0"

Then do:
abaqus job=XXX user=uel4nodeFbar.obj interactive

!*************************************************************************/


//#include <aba_for_c.h> since 6.13?
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cassert>
using namespace std;
#include "boost/multi_array.hpp"

//#include <Eigen/Dense>
//#include <Eigen/Eigenvalues>

typedef boost::multi_array<double, 4> Tensor4;
typedef boost::multi_array<double, 2> Tensor2;
typedef boost::multi_array<double, 3> MatrixD3;
typedef boost::multi_array<double, 2> MatrixD2;
typedef boost::multi_array<double, 1> VectorD;

// column-major for fortran array access
// h = row-height
#define FORT_ACCESS(A,r,c,h) ((A)[(r) + (c)*(h)])

void zeroM2(MatrixD2 & A);
void zeroT2(Tensor2 & A);
void zeroT4(Tensor4 & A);
void onem(Tensor2 & A);
double matlib_determinant(MatrixD2 & A);
double matlib_inverse(MatrixD2 & A, MatrixD2 & Ainv);
void matlib_transpose(MatrixD2 & A, MatrixD2 & AT);
void matlib_matmult(Tensor2 & A, Tensor2 & B, Tensor2 & C);

void UPE4(double * RHS, double * AMATRX, double * SVARS, double * ENERGY, int * NDOFEL, int * NRHS, int * NSVARS, double * PROPS, int * NPROPS, 
		  double * COORDS, int * MCRD, int * NNODE, double * U, double * DU, double * V, double * A, int * JTYPE, double * TIME, double * DTIME,
          int * KSTEP, int * KINC, int * JELEM, double * PARAMS, int * NDLOAD, int * JDLTYP, double * ADLMAG, double * PREDEF,
          int * NPREDF, int * LFLAGS,int * MLVARX,double * DDLMAG,int * MDLOAD,double * PNEWDT,int * JPROPS,int * NJPROP, double * PERIOD,
		  int nquad, int ndm, int nen);
		  
void compute_weight_upe4(int nIntPt, VectorD & w);
void compute_shape_upe4(int nquad, MatrixD2 & shp, MatrixD3 & dshp, VectorD & detj, MatrixD2 & x, bool shp_flag);

void elastic_constitutive_update(double * PROPS, int * NPROPS, Tensor2 & Fe_t, Tensor2 & F_t, Tensor2 & F_tau, Tensor2 & Fe_tau, Tensor2 & T_tau, Tensor4 & SpTanMod);
void neohookean_constitutive_update(double * PROPS,int * NPROPS, Tensor2 & Fe_t, Tensor2 & F_t, Tensor2 & F_tau, Tensor2 & Fe_tau, Tensor2 & T_tau, Tensor4 & SpTanMod);
	  
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
    if(LFLAGS[1]==1) {	
        // LFLAGS[1]=0 -> small disp.
        // LFLAGS[1]=1 -> large disp.
        cout << "Abaqus thinks you are doing" << endl;
        cout << "a small displacement analysis" << endl;
        cout << "go in and REMOVE nlgeom=yes" << endl;
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
		int nen   =  4;	

		UPE4(RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS, PROPS, NPROPS, 
		     COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, TIME, DTIME,
			 KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG, PREDEF,
			 NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT, JPROPS, NJPROP, PERIOD,
			 nquad, ndm, nen);
	
	} else {
		cout << "unrecognized element\n";
	}
      
}

void compute_weight_upe4(int nquad, VectorD & w)
{
	assert ( nquad == 1 || nquad == 4 );
	
	if ( nquad == 1 ) {
		
		w[0] = 2*2;
		
	} else if ( nquad == 4 ) {
	
		w[0] = 1.0*1.0;
        w[1] = 1.0*1.0;
		w[2] = 1.0*1.0; 
	    w[3] = 1.0*1.0;
		
	} else {
		cout << "compute_weight_upe4: shouldn't be here\n";	
	}
	
	
}

void compute_shape_upe4(int nquad, MatrixD2 & shp, MatrixD3 & dshp, VectorD & detj, MatrixD2 & x, bool shp_flag)
{
	
	assert ( nquad == 1 || nquad == 4 );
	
	int ndm  =  2;
	int nen  =  4;	
	
	MatrixD2 s(boost::extents[nquad][2]);
	MatrixD2 dnds(boost::extents[nen][ndm]); // temp
	MatrixD2 xjac(boost::extents[ndm][ndm]); // temp
	MatrixD2 dxds(boost::extents[ndm][ndm]); // temp
	double xsj;           		             // temp
	
	if ( nquad == 1 ) {
		
		s[0][0] = 0.0;
		s[0][1] = 0.0;	
	
	} else if ( nquad == 4 )  {
		
		double pt = 1./sqrt(3.);
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

	std::fill(  shp.origin(),  shp.origin() +  shp.num_elements(), 0.0 );
	std::fill( dshp.origin(), dshp.origin() + dshp.num_elements(), 0.0 );
	std::fill( detj.origin(), detj.origin() + detj.num_elements(), 0.0 );
	for (int lquad = 0; lquad < nquad; lquad++)
	{
		//cout << "gauss coord = (" << s[lquad][0] << "," << s[lquad][1] << endl;
		// ===================================
		// fill in shape functions N1,N2,N3
		// ===================================
		if ( shp_flag == true ) {
			shp[lquad][0] = 0.25*(1.0-s[lquad][0])*(1.0-s[lquad][1]); // N1 
			shp[lquad][1] = 0.25*(1.0+s[lquad][0])*(1.0-s[lquad][1]); // N2 
			shp[lquad][2] = 0.25*(1.0+s[lquad][0])*(1.0+s[lquad][1]); // N3 
			shp[lquad][3] = 0.25*(1.0-s[lquad][0])*(1.0+s[lquad][1]); // N4 		
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
		dnds[0][0] = -0.25*(1.0-s[lquad][1]); // dN1/dr
		dnds[0][1] = -0.25*(1.0-s[lquad][0]); // dN1/ds

		dnds[1][0] =  0.25*(1.0-s[lquad][1]); // dN2/dr
		dnds[1][1] = -0.25*(1.0+s[lquad][0]); // dN2/ds

		dnds[2][0] =  0.25*(1.0+s[lquad][1]); // dN3/dr
		dnds[2][1] =  0.25*(1.0+s[lquad][0]); // dN3/ds
		
		dnds[3][0] = -0.25*(1.0+s[lquad][1]); // dN4/dr
		dnds[3][1] =  0.25*(1.0-s[lquad][0]); // dN4/ds
		
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
		
		std::fill( dxds.origin(), dxds.origin() + dxds.num_elements(), 0.0 );
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
		
		assert( xsj > 1.e-12 );
		
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
		  int nquad, int ndm, int nen)
{
	int small_strain_flag = 0;
	int counter;  		  
	int nlSdv;  // number of state variables
	int ndf  = ndm;
	int nsdm = ndm*ndf;
	
	Tensor2 Fbar_t(boost::extents[3][3]);   double detFbar_t;     // previous converged Fbar
	Tensor2 Fbar_tau(boost::extents[3][3]); double detFbar_tau; // current Fbar
	
	Tensor2 F_t(boost::extents[3][3]);    double detF_t;       // previous converged def gradient
	Tensor2 Fe_t(boost::extents[3][3]);                     // previous converged elastic def gradient
	Tensor2 F_tau(boost::extents[3][3]);  double detF_tau;   // current def gradient
	Tensor2 Fe_tau(boost::extents[3][3]);                   // current elastic def gradient
	Tensor2 Fc_tau(boost::extents[3][3]); double detFc_tau; // current def gradient at centroid 'c' 
	Tensor2 Fc_t(boost::extents[3][3]);   double detFc_t;     // previous converged def gradient at centroid 'c' 
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
	
	// number of state variables
	nlSdv = *NSVARS; 
	
    // Initialize the residual and tangent matrices to zero.
    fill(Kuu.origin(),Kuu.origin()+Kuu.num_elements(),0.0);
    fill(Ru.origin(),Ru.origin()+Ru.num_elements(),0.0);
    *ENERGY=0.0; 
	
	// Body forces
    fill(body.origin(),body.origin()+body.num_elements(),0.0);
	
    // Obtain nodal displacements 
	counter = 0;
	for ( int i = 0; i < nen; ++i ) {
		for ( int j = 0; j < ndm; ++j ) {
			 u[i][j]   =  U[counter];
            du[i][j]   = DU[counter];
			uOld[i][j] = u[i][j]-du[i][j];
			//cout << "U[" << counter << "] = " << U[counter] << ", DU[" << counter << "] = " << DU[counter] << endl; 
			
			counter++;			
		}
	}
	
	// note: COORDS(dof,node)
	for ( int i = 0; i < nen; ++i ) {
		for ( int j = 0; j < ndm; ++j ) {
			 coords[i][j] = FORT_ACCESS(COORDS,j,i,ndm);
			coordsC[i][j] = coords[i][j];// + small_strain_flag*u[i][j];		
			//cout << "u[" << i << "][" << j << "] = " << u[i][j] 
			 //    << ", coords[" << i << "][" << j << "] = " << coords[i][j]
				// << ", coordsC[" << i << "][" << j << "] = " << coordsC[i][j] << endl;
			
				
		}
	}

	// ----------------------------------------------------------------
	// compute shape functions and derivatives
    if(nen == 4) {
		
		compute_shape_upe4(    1, shp0,  dshp0,  detj0,  coords, true ); // at centroid (0), reference, for F-bar method
		compute_shape_upe4(    1, shp0, dshpC0, detjC0, coordsC, false); // at centroid (0), current (C), for F-bar method
        compute_shape_upe4(nquad,  shp,   dshp,   detj,  coords, true ); // at all gauss points, reference
        compute_shape_upe4(nquad,  shp,  dshpC,  detjC, coordsC, false); // at all gauss points, current (C)
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
	for ( int i=0; i < ndm; ++i ) {
		for ( int j=0; j < ndm; ++j ) {
			for ( int k=0; k < nen; ++k ) {
				Fc_tau[i][j] += u[k][i]*dshp0[0][k][j];
				Fc_t[i][j] += uOld[k][i]*dshp0[0][k][j];
			}
		}
	}
	// modify for plane-strain
	Fc_tau[2][2] = 1.0;
	Fc_t[2][2] = 1.0;
	
	// 2D plane-strain implementation detF
	detFc_t   = Fc_t[0][0]*Fc_t[1][1] - Fc_t[0][1]*Fc_t[1][0];
	detFc_tau = Fc_tau[0][0]*Fc_tau[1][1] - Fc_tau[0][1]*Fc_tau[1][0];
	//
	// With the deformation gradient known at the element centriod
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
			onem(F_t);
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
	
        // Modify the deformation gradient for the `F-bar' method
        // only when using the 4 node fully integrated linear
        // element, do not use the `F-bar' method for any other element
		//
		// Note: LFLAGS[1]==1 means large-displacement analysis
        //
		/*if ((nen==4)&&(nquad==4)&&(LFLAGS[1]==1)&&(small_strain_flag==1)) 
		{
			//  2D plane-strain implementation
            detF_t   = F_t[0][0]*F_t[1][1] - F_t[0][1]*F_t[1][0];
			detF_tau = F_tau[0][0]*F_tau[1][1] - F_tau[0][1]*F_tau[1][0];
			double fac_tau = pow((detFc_tau/detF_tau),0.5);
			double fac_t = pow((detFc_t/detF_t),0.5);
			
			for ( int i = 0; i < 3; ++i ) {
				for ( int j = 0; j < 3; ++j ) { 
					F_tau[i][j] = fac_tau*F_tau[i][j];
					F_t[i][j]   = fac_t*F_t[i][j];
				}
			}	
        }*/
		
		// Perform the constitutive time integration at this integ. point 
		// Input:
		//   PROPS  = material properties
		//   NPROPS = number of properties
		//   Fe_t   = previously converged elastic def gradient (Fe-bar previous)
		//   F_t    = previously converged total def gradient (F-bar previous)
		//   F_tau  = current total def gradient (F-bar current)
		//
		// Output:
		//   Fe_tau   = current elastic def gradient (Fe-bar current)
		//   T_tau    = current Cauchy stress
		//   SpTanMod = current tangent moduli
		//
		//neohookean_constitutive_update(PROPS,NPROPS,Fe_t,F_t,F_tau,Fe_tau,T_tau,SpTanMod);
		elastic_constitutive_update(PROPS,NPROPS,Fe_t,F_t,F_tau,Fe_tau,T_tau,SpTanMod);

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
		
		/*cout << "Amat = " << endl;
		for ( int i = 0; i < nsdm; ++i ) {
			for ( int j = 0; j < nsdm; ++j ) {
				cout.width(10); 
				cout << std::scientific << std::right << Amat[i][j] << " ";
			}
			cout << endl;
		}*/
			
        zeroM2(Qmat); // tangent moduli at current integration point
		Qmat[0][0] = 0.5*(Amat[0][0]+Amat[0][3]) - 0.5*T_tau[0][0];
		Qmat[1][0] = 0.5*(Amat[1][0]+Amat[1][3]) - 0.5*T_tau[0][1];
		Qmat[2][0] = 0.5*(Amat[2][0]+Amat[2][3]) - 0.5*T_tau[0][1];
		Qmat[3][0] = 0.5*(Amat[3][0]+Amat[3][3]) - 0.5*T_tau[1][1];
		Qmat[0][3] = 0.5*(Amat[0][0]+Amat[0][3]) - 0.5*T_tau[0][0];
		Qmat[1][3] = 0.5*(Amat[1][0]+Amat[1][3]) - 0.5*T_tau[0][1];
		Qmat[2][3] = 0.5*(Amat[2][0]+Amat[2][3]) - 0.5*T_tau[0][1];
		Qmat[3][3] = 0.5*(Amat[3][0]+Amat[3][3]) - 0.5*T_tau[1][1];


		// Swansea style
		for (int a=0; a<nen; a++) { // loop over nodes
			for(int i=0; i<ndf; i++){ // loop over dofs of node
				for(int b=0; b<nen; b++) { // loop over nodes
					for(int j=0; j<ndf; j++) { // loop over dofs of node					
						for(int I=0; I<nsdm; I++) { // loop over nsdm == ndm*ndf (tangent)
							for(int J=0; J<nsdm; J++) { // loop over nsdm == ndm*ndf (tangent)
								Kuu[a*ndf+i][b*ndf+j]+=(  Gmat[I][ndf*a+i]*Amat[I][J]*Gmat[J][ndf*b+j]  )*detjC[lquad]*w[lquad];
														//+ Gmat[I][ndf*a+i]*Qmat[I][J]*( G0mat[J][ndf*b+j]-Gmat[J][ndf*b+j] )  )*detjC[lquad]*w[lquad];
							}
						}
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


	// Assemble the element level residual
	for ( int a = 0; a < nen; a++ ) {
		for(int i = 0; i < ndf; i++ ) { // loop over dofs of node
			RHS[a*ndf+i] = Ru[a][i];
			//FORT_ACCESS(RHS,i,a,ndf) = Ru[a][i]; 
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
	// amatrx[row + col*8] = Kuu[row][col]
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

	cout << "ke = " << endl;
	for ( int i = 0; i < nen*ndf; ++i ) {
		for ( int j = 0; j < nen*ndf; ++j ) {
			cout.width(20);
			cout.precision(12);
			cout << std::right << Kuu[i][j] << " ";
			//FORT_ACCESS(AMATRX,i,j,nen*ndf)=Kuu[i][j];
		}
		cout << endl;
	}
	//cin.get();
		
		
}

void elastic_constitutive_update(double * PROPS, int * NPROPS, Tensor2 & Fe_t, Tensor2 & F_t, Tensor2 & F_tau, Tensor2 & Fe_tau, Tensor2 & T_tau, Tensor4 & SpTanMod)
{
	Tensor2 Iden(boost::extents[3][3]);
	Tensor2 Eps(boost::extents[3][3]);
	double E,nu,theta,fac1,fac2;
	
	// identity tensor
    onem(Iden);
	
    // Obtain material properties
    E  = PROPS[0];
    nu = PROPS[1];
	
	// Small strain tensor  
	zeroT2(Eps);
	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
			Eps[i][j] = F_tau[i][j] - Iden[i][j];
	
	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
			Eps[i][j] = 0.5*(Eps[i][j]+Eps[j][i]);
	
	//cout << "Eps[2][2] = " << Eps[2][2] << endl;
	
	theta = 0;
	for ( int i = 0; i < 3; ++i )
		theta += Eps[i][i];
	
	// compute Cauchy stress (see lecture notes of Ortiz)
	zeroT2(T_tau);
	fac1 = E/(1.0+nu);
    fac2 = theta*nu/(1.0-2.0*nu);
    //T_tau = fac1*( eps + fac2*Iden );
	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
			T_tau[i][j] = fac1*( Eps[i][j] + fac2*Iden[i][j]  );
	
	  
	// compute spatial tangent moduli (see lecture notes of Ortiz)
	zeroT4(SpTanMod);	
	fac1 = E/(2.0*(1.0+nu));
    fac2 = E*nu/((1.0+nu)*(1.0-2.0*nu));
    //SpUUMod(i,j,k,l) = fac1*(Iden(i,l)*Iden(j,k)+Iden(i,k)*Iden(j,l)) + fac2*Iden(i,j)*Iden(k,l)  	
	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
			for ( int k = 0; k < 3; ++k )
				for ( int l = 0; l < 3; ++l )
					SpTanMod[i][j][k][l] = fac1*( Iden[i][l]*Iden[j][k] + Iden[i][k]*Iden[j][l]  ) + fac2*Iden[i][j]*Iden[k][l];
	
}

/*void neohookean_constitutive_update(double * PROPS, int * NPROPS, Tensor2 & Fe_t, Tensor2 & F_t, Tensor2 & F_tau, Tensor2 & Fe_tau, Tensor2 & T_tau, Tensor4 & SpTanMod)
{
	// hyperelastic compressible neo-hookean
	Tensor2 Iden(boost::extents[3][3]);
	Tensor2 F_tauT(boost::extents[3][3]);
	Tensor2 B(boost::extents[3][3]);
	//Tensor2 Cinv(boost::extents[3][3]);
	//Tensor2 S(boost::extents[3][3]);
	//Tensor2 P(boost::extents[3][3]);
	//Tensor4 LagTanMod(boost::extents[3][3][3][3]);
	double lambda,mu;
	double J,logJ,detC;
	//double coef;
	
	// identity tensor
    onem(Iden);
	
    // Obtain material properties
    lambda = PROPS[0];
    mu     = PROPS[1];
	cout << "lambda = " << lambda << endl;
	cout << "mu = " << mu << endl;
	
	J    = matlib_determinant(F_tau);
	logJ = log(J);
	
	// transpose of F
	matlib_transpose(F_tau,F_tauT);
	
	// compute left Cauchy-Green tensor B=F*F^T
	matlib_matmult(F_tau, F_tauT, B);
	
	// compute Cauchy stress (see lecture notes of Ortiz)
	zeroT2(T_tau);
	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
			T_tau[i][j] += (mu*B[i][j]+(lambda*logJ-mu)*Iden[i][j])/J;
	
	
	// compute spatial tangent moduli (see lecture notes of Ortiz)
	zeroT4(SpTanMod);				
	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
			for ( int k = 0; k < 3; ++k )
				for ( int l = 0; l < 3; ++l )
					SpTanMod[i][j][k][l] += (  lambda*Iden[i][j]*Iden[k][l] + (mu-lambda*logJ)*(Iden[i][k]*Iden[j][l]+Iden[i][l]*Iden[j][k])  )/J + Iden[i][k]*T_tau[j][l];
		
	
	// for neo-hookean total def gradient is elastic
	for ( int i = 0; i < 3; ++i ) {
		for ( int J = 0; J < 3; ++J ) {
			Fe_tau[i][J] = F_tau[i][J];
		}
	}
	
}*/

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

