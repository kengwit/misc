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
!  *User Element,Nodes=4,Type=U1,Iproperties=1,Properties=2,Coordinates=2,Variables=36,Unsymm
!  1,2
!
! Note: 9 state variables are used in this element to store previously 
! converged elastic deformation gradient (Fe-bar)
!
! No integer properties (jprops) used for this element - we only use nInt = 4 (full integration) 
!
! Material Properties Vector
! --------------------------------------------------------------
! lambda  = props[0]  ! lambda
! mu      = props[1]  ! mu
! nu      = lambda/(2*(lambda+mu))
! nlSdv   = jprops[0] ! number of state variables PER INTEGRATION POINT 

The C++ umat must be compiled externally using a compiler supported by ABAQUS, producing an object file .obj. The compiler command is

cl /EHsc /c /Fo uel4nodeFbarBlitz.cpp /I"C:\Users\KWL\Desktop\Research\libraries\blitz-0.10"

Then do:
abaqus job=XXX user=uel4nodeFbar.obj interactive

!*************************************************************************/


//#include <aba_for_c.h> since 6.13?
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cassert>
#include <cmath>
using namespace std;

#include <blitz/array.h>
typedef blitz::Array<double,4> MatD4;
typedef blitz::Array<double,3> MatD3;
typedef blitz::Array<double,2> MatD2;
typedef blitz::Array<double,1> VecD;

// column-major for fortran array access
// h = row-height
#define FORT_ACCESS(A,r,c,h) ((A)[(r) + (c)*(h)])

void onem(MatD2 A)
{
	A = 0.0;
	for ( int i = 0; i < 3; i++ ) A(i,i) = 1.0;
}

void UPE4(double * RHS, double * AMATRX, double * SVARS, double * ENERGY, int * NDOFEL, int * NRHS, int * NSVARS, double * PROPS, int * NPROPS, 
		  double * COORDS, int * MCRD, int * NNODE, double * U, double * DU, double * V, double * A, int * JTYPE, double * TIME, double * DTIME,
          int * KSTEP, int * KINC, int * JELEM, double * PARAMS, int * NDLOAD, int * JDLTYP, double * ADLMAG, double * PREDEF,
          int * NPREDF, int * LFLAGS,int * MLVARX,double * DDLMAG,int * MDLOAD,double * PNEWDT,int * JPROPS,int * NJPROP, double * PERIOD,
		  int nquad, int ndm);
		  
void neohookean_constitutive_update(double * PROPS, int * NPROPS, MatD2 Fe_t, MatD2 F_t, MatD2 F_tau, MatD2 Fe_tau, MatD2 T_tau, MatD4 SpTanMod);
	  
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


void get_xintpt(MatD2 xi, VecD w, int nintpts)
{
	// Initialize
	w  = 0.0;
	xi = 0.0;

	if ( nintpts == 1 ) {
			
		w(0) = 4.0;
      
		xi(0,0) = 0.0;
		xi(0,1) = 0.0;
		
	} else if ( nintpts == 4 ) {
      
		// Gauss weights
		//
		w(0) = 1.0;
		w(1) = 1.0;
		w(2) = 1.0;
		w(3) = 1.0;

		// Gauss pt locations in master element
		xi(0,0) = -sqrt(1.0/3.0);
		xi(0,1) = -sqrt(1.0/3.0);

		xi(1,0) = sqrt(1.0/3.0);
		xi(1,1) = -sqrt(1.0/3.0);

		xi(2,0) = sqrt(1.0/3.0);
		xi(2,1) = sqrt(1.0/3.0);

		xi(3,0) = -sqrt(1.0/3.0);
		xi(3,1) = sqrt(1.0/3.0);

    } else {
		cout << "wrong number of integration points\n";
		exit(1);
	}
  	
}
void calcShape2DLinear(int nintpts, MatD2 xi_int, int intpt, VecD sh, MatD2 dshxi )
{
	double xi  = xi_int(intpt,0);
    double eta = xi_int(intpt,1);
    double one = 1.0;
	double fourth = 0.25;
      
    // The shape functions
	sh(0) = fourth*(one - xi)*(one - eta);
	sh(1) = fourth*(one + xi)*(one - eta);
	sh(2) = fourth*(one + xi)*(one + eta);
	sh(3) = fourth*(one - xi)*(one + eta);

      
    // The first derivatives    
	dshxi(0,0) = -fourth*(one - eta);
	dshxi(0,1) = -fourth*(one - xi);

	dshxi(1,0) = fourth*(one - eta);
	dshxi(1,1) = -fourth*(one + xi);

	dshxi(2,0) = fourth*(one + eta);
	dshxi(2,1) = fourth*(one + xi);

	dshxi(3,0) = -fourth*(one + eta);
	dshxi(3,1) = fourth*(one - xi);


}

void mapShape2Da(int nNode, MatD2 dshxi, MatD2 coords, MatD2 dsh, double & detMapJ)
{
	MatD2 mapJ(2,2);
	MatD2 mapJ_inv(2,2);
	
    // Calculate the mapping Jacobian matrix:
	mapJ = 0.0;
	for ( int i = 0; i < 2; i++ )
		for ( int j = 0; j < 2; j++ )
			for ( int k = 0; k < nNode; k++ )
				mapJ(i,j) += dshxi(k,i)*coords(k,j);
	  

	// ==============================
	// Compute determinant of Jacobian matrix
	// ==============================
	detMapJ = mapJ(0,0)*mapJ(1,1)-mapJ(0,1)*mapJ(1,0);
	
	if ( detMapJ < 1.e-20 )
	{
		cout << "detMapJ = " << detMapJ << endl;
		cout << "determinant Jacobian close to zero or negative. Exiting ...\n";
		exit(1);
	}
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
	mapJ_inv(0,0) = mapJ(1,1)/detMapJ; // dr/dx
	mapJ_inv(0,1) =-mapJ(0,1)/detMapJ; // ds/dx
	mapJ_inv(1,0) =-mapJ(1,0)/detMapJ; // dr/dy
	mapJ_inv(1,1) = mapJ(0,0)/detMapJ; // ds/dy

	// ==============================
	// compute and store shape function 
	// derivatives wrt to real coordinates
	// ==============================
	// just think of Jac^inv * {dNi/dr}
	//                         {dNi/ds} 	
	dsh = 0.0;
	for ( int a = 0; a < nNode; ++a )
		for ( int i = 0; i < 2; ++i )
			for ( int j = 0; j < 2; ++j )
				dsh(a,i) += mapJ_inv(i,j)*dshxi(a,j);



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
	int nlSdv = JPROPS[0];  // number of state variables PER INTEGRATION POINT
	
	MatD2 Fbar_t(3,3);   double detFbar_t;     // previous converged Fbar
	MatD2 Fbar_tau(3,3); double detFbar_tau; // current Fbar
	
	MatD2 F_t(3,3);      double detF_t;       // previous converged def gradient
	MatD2 Fe_t(3,3);                     // previous converged elastic def gradient
	
	MatD2 F_tau(3,3);    double detF_tau;   // current def gradient
	MatD2 Fe_tau(3,3);                   // current elastic def gradient
	
	MatD2 Fc_tau(3,3);   double detFc_tau; // current def gradient at centroid 'c' 
	MatD2 Fc_t(3,3);     double detFc_t;     // previous converged def gradient at centroid 'c' 
	
	MatD2 T_tau(3,3);
	
	MatD4 SpTanMod(3,3,3,3); 
	
	// ------------------------------------------------------
	// shape functions & derivatives at centroid
	// ------------------------------------------------------
	MatD2 xi0(1,2);  
	
	VecD sh0(nen);         
	MatD2 dshxi0(nen,ndm);
	
	MatD2 dsh0(nen,ndm);   // ref config, used to calculate deformation gradient at centroid
	double  detMapJ0;      // ref config, not used
	
	MatD2 dshC0(nen,ndm);  // current config, used to calculate G0mat
	double  detMapJ0C;     // current config, not used
	
	VecD w0(1);            // not used
	
	// ------------------------------------------------------
	// shape functions & derivatives at normal gauss points
	// ------------------------------------------------------
	MatD2 xi(nquad,2);
	MatD2 dshxi(nen,ndm);
	
	VecD sh(nen);        // shape functions
	
	MatD2 dsh(nen,ndm);  // ref config, used to calculate deformation gradient
	double detMapJ;      // ref config, not used
		
	MatD2 dshC(nen,ndm); // current config, used for Bmat, Gmat
	double detMapJC;     // current config, used for Bmat, Gmat
	
	VecD  w(nquad);      // weights
	  	
	
	// ------------------------------------------------------
	// coordinates & displacements
	// ------------------------------------------------------
	MatD2 u(nen,ndm);
	MatD2 du(nen,ndm);
	MatD2 uOld(nen,ndm);
	MatD2 coords(nen,ndm);
	MatD2 coordsC(nen,ndm);
	
	// ------------------------------------------------------
	// Voigt style quantities
	// ------------------------------------------------------
	MatD2 Kuu(nen*ndm,nen*ndm);
	MatD2 Ru(nen,ndm);
	VecD  body(ndm);
	
	VecD  Smat(3);
	MatD2 Bmat(3,ndm*nen);
	MatD2 Gmat(4,ndm*nen);
	MatD2 G0mat(4,ndm*nen);
	MatD2 Qmat(4,4);
	MatD2 Amat(4,4);
		
	// Initialize the residual and tangent matrices to zero.
    Kuu = 0.0;
    Ru = 0.0;
    *ENERGY=0.0; 
	
	// Body forces
	body = 0.0;
	
	// Obtain nodal displacements 
	counter = 0;
	for ( int a = 0; a < nen; ++a ) {
		for ( int j = 0; j < ndm; ++j ) {
			 u(a,j)   =  U[counter];
            du(a,j)   = DU[counter];
			uOld(a,j) = u(a,j)-du(a,j);
			counter++;			
		}
	}
	
	// COORDS is column major
	// 0 2 4
	// 1 3 5
	for ( int a = 0; a < nen; ++a ) {
		for ( int j = 0; j < ndm; ++j ) {
			 coords(a,j) = COORDS[j+a*ndm];
			coordsC(a,j) = coords(a,j) + u(a,j);		
			
		}
	}
	
	
	/*----------------------------------------------------------------
	! 
	! Take this opportunity to perform calculations at the element
	!  centroid.  Get the deformation gradient for use in the
	!  `F-bar' method.
	!
	! Reference for the F-bar method:
	!  de Souza Neto, E.A., Peric, D., Dutko, M., Owen, D.R.J., 1996.
	!  Design of simple low order finite elements for large strain
	!  analysis of nearly incompressible solids. International Journal
	!  of Solids and Structures, 33, 3277-3296.
	!
	!
	! Obtain shape functions and their local gradients at the element
	!  centriod, that means xi=eta=zeta=0.0, and nIntPt=1
	*/
	get_xintpt(xi0,w0,1);
	calcShape2DLinear(1,xi0,0,sh0,dshxi0);
	mapShape2Da(nen,dshxi0,coords,dsh0,detMapJ0);
	mapShape2Da(nen,dshxi0,coordsC,dshC0,detMapJ0C);

	// ==================================================================
	// Calculate the deformation gradient at the element centriod
    // at the the begining and end of the increment for use in 
    // the `F-bar' method. `Tau' represents the end of the increment
    // and `t' the previous increment.
    onem(Fc_tau);
    onem(Fc_t);
	
	for ( int i=0; i < ndm; ++i ) {
		for ( int j=0; j < ndm; ++j ) {
			for ( int a=0; a < nen; ++a ) {
				Fc_tau(i,j) +=    u(a,i)*dsh0(a,j);
				Fc_t(i,j)   += uOld(a,i)*dsh0(a,j);
			}
		}
	}
	
	
	// modify for plane-strain
	Fc_tau(2,2) = 1.0;
	Fc_t(2,2) = 1.0;
	
	// 2D plane-strain implementation detF
	detFc_t   =     Fc_t(0,0)*Fc_t(1,1) - Fc_t(0,1)*Fc_t(1,0);
	detFc_tau = Fc_tau(0,0)*Fc_tau(1,1) - Fc_tau(0,1)*Fc_tau(1,0);
	//
	// With the deformation gradient known at the element centroid
	// we are now able to implement the `F-bar' method later
	// ==================================================================
	
	int jj;
	double * svars  = SVARS; // current block of state variables
	// loop through integration points
	get_xintpt(xi,w,nquad);
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
					Fe_t(i,j) = svars[jj];
					jj++;
				}
			}
			
        }

        calcShape2DLinear(nquad,xi,lquad,sh,dshxi);
        mapShape2Da(nen,dshxi,coords,dsh,detMapJ);
        mapShape2Da(nen,dshxi,coordsC,dshC,detMapJC);

		// ==================================================================
		// Obtain, and modify the deformation gradient at this integration
        //  point.  Modify the deformation gradienet for use in the `F-bar'
        //  method.  Also, take care of plane-strain.
		//
        onem(F_tau); 
        onem(F_t);   
	
		for ( int i=0; i < ndm; ++i ) {
			for ( int j=0; j < ndm; ++j ) {
				for ( int a=0; a < nen; ++a ) {
					F_tau(i,j) +=    u(a,i)*dsh(a,j);
					F_t(i,j)   += uOld(a,i)*dsh(a,j);
				}
			}
		}
		
		
		// modify for plane-strain
		F_tau(2,2) = 1.0;
		F_t(2,2) = 1.0;
		
		
        // Modify the deformation gradient for the `F-bar' method
        // only when using the 4 node fully integrated linear
        // element, do not use the `F-bar' method for any other element
		//
		// Note: LFLAGS[1]==1 means large-displacement analysis
        //
		if ((nen==4)&&(nquad==4)&&(LFLAGS[1]==1)) 
		{
			// 2D plane-strain implementation detF
			detF_t   =     F_t(0,0)*F_t(1,1) - F_t(0,1)*F_t(1,0);
			detF_tau = F_tau(0,0)*F_tau(1,1) - F_tau(0,1)*F_tau(1,0);
		
			double fac_t =  pow((detFc_t/detF_t),0.5);
			double fac_tau = pow((detFc_tau/detF_tau),0.5);
			
			// note below the factor does not multiply F33 (which remains at 1.0)
			for ( int i = 0; i < ndm; ++i ) {
				for ( int j = 0; j < ndm; ++j ) { 
					F_t(i,j)   = fac_t*F_t(i,j);
					F_tau(i,j) = fac_tau*F_tau(i,j);					
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
		//   Fe_t   = previously converged elastic def gradient (Fe-bar previous)
		//   F_t    = previously converged total def gradient (F-bar previous)
		//   F_tau  = current total def gradient (F-bar current)
		//
		// Output:
		//   Fe_tau   = current elastic def gradient (Fe-bar current)
		//   T_tau    = current Cauchy stress
		//   SpTanMod = current tangent moduli
		//
		neohookean_constitutive_update(PROPS,NPROPS,Fe_t,F_t,F_tau,Fe_tau,T_tau,SpTanMod);

		// Compute/update the displacement residual vector
        Smat(0) = T_tau(0,0);
		Smat(1) = T_tau(1,1);
		Smat(2) = T_tau(0,1);
		
		Bmat = 0.0;
		for ( int a=0; a < nen; a++ )
		{
			Bmat(0,ndm*a  ) = dshC(a,0); // dNi/dx
			Bmat(1,ndm*a+1) = dshC(a,1); // dNi/dy
			Bmat(2,ndm*a  ) = dshC(a,1); // dNi/dy
			Bmat(2,ndm*a+1) = dshC(a,0); // dNi/dx
		}
		
		// Swansea style
		for ( int a = 0; a < nen; a++ ) { // loop over nodes
			for ( int i = 0; i < ndm; ++i ) { // loop over degrees of freedom
				for ( int j = 0; j < 3; ++j ) { // loop over stress terms
					Ru(a,i) -= Bmat(j,ndm*a+i)*Smat(j)*detMapJC*w(lquad);
				}
			} 
		} // end loop over nodes 		

		// Compute/update the displacement tangent matrix
        
        Gmat = 0.0;
		for ( int a=0; a < nen; a++ )
		{
			Gmat(0,ndm*a  ) = dshC(a,0); // dNi/dx
			Gmat(1,ndm*a+1) = dshC(a,0); // dNi/dx
			Gmat(2,ndm*a  ) = dshC(a,1); // dNi/dy
			Gmat(3,ndm*a+1) = dshC(a,1); // dNi/dy
		}
		
		G0mat = 0.0; // tangent moduli at centroid (note: only at centroid, index = 0)
		for ( int a=0; a < nen; a++ )
		{
			G0mat(0,ndm*a  ) = dshC0(a,0); // dNi/dx at centroid
			G0mat(1,ndm*a+1) = dshC0(a,0); // dNi/dx at centroid
			G0mat(2,ndm*a  ) = dshC0(a,1); // dNi/dy at centroid
			G0mat(3,ndm*a+1) = dshC0(a,1); // dNi/dy at centroid
		}
			
		Amat = 0.0; // tangent moduli at current integration point
		
		Amat(0,0) = SpTanMod(0,0,0,0); // 1111 -> 11 -> node1, dNi/dx (must jibe with Gmat above)
		Amat(0,1) = SpTanMod(0,0,1,0); // 1121 -> 21 -> node2, dNi/dx
		Amat(0,2) = SpTanMod(0,0,0,1); // 1112 -> 12 -> node1, dNi/dy
		Amat(0,3) = SpTanMod(0,0,1,1); // 1122 -> 22 -> node2, dNi/dy
		
		Amat(1,0) = SpTanMod(1,0,0,0); // 2111 -> 11 -> node1, dNi/dx (must jibe with Gmat above)
		Amat(1,1) = SpTanMod(1,0,1,0); // 2121 -> 21 -> node2, dNi/dx
		Amat(1,2) = SpTanMod(1,0,0,1); // 2112 -> 12 -> node1, dNi/dy
		Amat(1,3) = SpTanMod(1,0,1,1); // 2122 -> 22 -> node2, dNi/dy
		
		Amat(2,0) = SpTanMod(0,1,0,0); // 1211 -> 11 -> node1, dNi/dx (must jibe with Gmat above)
		Amat(2,1) = SpTanMod(0,1,1,0); // 1221 -> 21 -> node2, dNi/dx
		Amat(2,2) = SpTanMod(0,1,0,1); // 1212 -> 12 -> node1, dNi/dy
		Amat(2,3) = SpTanMod(0,1,1,1); // 1222 -> 22 -> node2, dNi/dy
		
		Amat(3,0) = SpTanMod(1,1,0,0); // 2211 -> 11 -> node1, dNi/dx (must jibe with Gmat above)
		Amat(3,1) = SpTanMod(1,1,1,0); // 2221 -> 21 -> node2, dNi/dx
		Amat(3,2) = SpTanMod(1,1,0,1); // 2212 -> 12 -> node1, dNi/dy
		Amat(3,3) = SpTanMod(1,1,1,1); // 2222 -> 22 -> node2, dNi/dy
		
        Qmat = 0.0; // tangent moduli at current integration point
		
		Qmat(0,0) = 0.5*( Amat(0,0)+Amat(0,3) ) - 0.5*T_tau(0,0);
		Qmat(1,0) = 0.5*( Amat(1,0)+Amat(1,3) ) - 0.5*T_tau(0,1);
		Qmat(2,0) = 0.5*( Amat(2,0)+Amat(2,3) ) - 0.5*T_tau(0,1);
		Qmat(3,0) = 0.5*( Amat(3,0)+Amat(3,3) ) - 0.5*T_tau(1,1);

		Qmat(0,3) = 0.5*( Amat(0,0)+Amat(0,3) ) - 0.5*T_tau(0,0);
		Qmat(1,3) = 0.5*( Amat(1,0)+Amat(1,3) ) - 0.5*T_tau(0,1);
		Qmat(2,3) = 0.5*( Amat(2,0)+Amat(2,3) ) - 0.5*T_tau(0,1);
		Qmat(3,3) = 0.5*( Amat(3,0)+Amat(3,3) ) - 0.5*T_tau(1,1);
	
		for(int i=0; i < nen*ndf; i++){ 
			for(int j=0; j < nen*ndf; j++) { 					
				for(int I=0; I < nsdm; I++) { 
					for(int J=0; J < nsdm; J++) { 
						Kuu(i,j) += detMapJC*w(lquad)*Gmat(I,i)*( 
									Amat(I,J)*Gmat(J,j)
								+   Qmat(I,J)*(G0mat(J,j)-Gmat(J,j))
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
				svars[jj] = Fe_tau(i,j);
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
			RHS[a*ndf+i] = Ru(a,i);
		}
	}
	
	// Assemble the element level tangent matrix
	int row,col;
	for (int a=0; a<nen; a++) { // loop over nodes
		for(int i=0; i<ndf; i++){ // loop over dofs of node
			for(int b=0; b<nen; b++) { // loop over nodes
				for(int j=0; j<ndf; j++) { // loop over dofs of node					
					row = a*ndf+i;
					col = b*ndf+j;
					FORT_ACCESS(AMATRX,row,col,nen*ndf)=Kuu(row,col);
				}
			}
		}
	}
				
}

void neohookean_constitutive_update(double * PROPS, int * NPROPS, MatD2 Fe_t, MatD2 F_t, MatD2 F_tau, MatD2 Fe_tau, MatD2 T_tau, MatD4 SpTanMod)
{
	// hyperelastic compressible neo-hookean (see Ortiz lecture notes)
	// strain energy function PER UNIT UNDEFORMED VOLUME = W = 0.5*mu ( tr(C) - 3 ) + f(J)
	// where f(J) = 0.5*lambda*(logJ)^2 - mu*logJ
	// --> f'(J) = (lambda*logJ-mu)/J
	// --> f''(J) = (lambda-lambda*logJ+mu)/J^2
	blitz::firstIndex  i;
	blitz::secondIndex j;
	blitz::thirdIndex  k;
	blitz::fourthIndex l;
	
	MatD2 Iden(3,3);
	MatD2 B(3,3);
	double lambda,mu;
	double detF;
	
	// identity tensor
    onem(Iden);
	
    // Obtain material properties
    lambda = PROPS[0];
    mu     = PROPS[1];
	
	detF = F_tau(0,0)*F_tau(1,1) - F_tau(0,1)*F_tau(1,0);
	
	assert(detF > 1.e-20);
	
	// compute left Cauchy-Green tensor B=F*F^T
	B = blitz::sum(F_tau(i,k)*F_tau(j,k),k);
	/*B=0.0;
	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
			for ( int k = 0; k < 3; ++k )
				B(i,j) += F_tau(i,k)*F_tau(j,k);*/
		
	// compute Cauchy stress (see lecture notes of Ortiz)
	// eq. 4.7.22c
	T_tau = ( mu*B(i,j) + (lambda*log(detF)-mu)*Iden(i,j) )/detF;
	//T_tau = ( mu*B + (lambda*log(detF)-mu)*Iden )/detF;
	
	//T_tau = 0.0;
	//for ( int i = 0; i < 3; ++i )
	//	for ( int j = 0; j < 3; ++j )
	//		T_tau(i,j) = ( mu*B(i,j) + (lambda*log(detF)-mu)*Iden(i,j) )/detF;
	
	//cout << "B in umat " << endl << B << endl;
	//cout << "Tau in umat " << endl << T_tau << endl;
	
	// compute spatial tangent moduli (see lecture notes of Ortiz)
	// eq. 4.7.26b
	//SpTanMod = 0.0;
	//for ( int i = 0; i < 3; ++i )
	//	for ( int j = 0; j < 3; ++j )
	//		for ( int k = 0; k < 3; ++k )
	//			for ( int l = 0; l < 3; ++l )
	//				SpTanMod(i,j,k,l) = (  lambda*Iden(i,j)*Iden(k,l) + (mu-lambda*log(detF))*(Iden(i,k)*Iden(j,l)+Iden(i,l)*Iden(j,k))  )/detF + Iden(i,k)*T_tau(j,l);
	
	MatD4 t1(3,3,3,3);
	MatD4 t2(3,3,3,3);
	MatD4 t3(3,3,3,3);
	MatD4 t4(3,3,3,3);
	
	t1 = lambda*Iden(i,j)*Iden(k,l);
	t2 = (mu-lambda*log(detF))*( Iden(i,k)*Iden(j,l) + Iden(i,l)*Iden(j,k) );
	t3 = Iden(i,k)*T_tau(j,l);
	SpTanMod = ( t1 + t2 )/detF + t3;
	
	// for neo-hookean total def gradient is elastic
	Fe_tau = F_tau;
	
	//for ( int i = 0; i < 3; ++i ) {
	//	for ( int j = 0; j < 3; ++j ) {
	//		Fe_tau(i,j) = F_tau(i,j);
	//	}
	//}
	
}
