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
! lambda  = props(1)  ! lambda
! mu      = props(2)  ! mu
! nu      = lambda/(2*(lambda+mu))

The C++ umat must be compiled externally using a compiler supported by ABAQUS, producing an object file .obj. The compiler command is

cl /EHsc /c /Fo uel4nodeFbarEPBlitzIncomplete.cpp /I"C:\Users\KWL\Desktop\Research\libraries\eigen-3.2.4" /I"C:\Users\KWL\Desktop\Research\libraries\blitz-0.10"

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

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
typedef Eigen::MatrixXd EgMatD2;
typedef Eigen::VectorXd EgVecD;

#include <blitz/array.h>
typedef blitz::Array<double,4> MatD4;
typedef blitz::Array<double,3> MatD3;
typedef blitz::Array<double,2> MatD2;
typedef blitz::Array<double,1> VecD;

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
void solve_eigenvalue_3x3(MatD2 AMat, VecD BePr, MatD2 Nvec);

double matlib_determinant(MatD2 A);
double matlib_inverse(MatD2 A, MatD2 Ainv);
void onem(MatD2 A);

void UPE4(double * RHS, double * AMATRX, double * SVARS, double * ENERGY, int * NDOFEL, int * NRHS, int * NSVARS, double * PROPS, int * NPROPS, 
		  double * COORDS, int * MCRD, int * NNODE, double * U, double * DU, double * V, double * A, int * JTYPE, double * TIME, double * DTIME,
          int * KSTEP, int * KINC, int * JELEM, double * PARAMS, int * NDLOAD, int * JDLTYP, double * ADLMAG, double * PREDEF,
          int * NPREDF, int * LFLAGS,int * MLVARX,double * DDLMAG,int * MDLOAD,double * PNEWDT,int * JPROPS,int * NJPROP, double * PERIOD,
		  int nquad, int ndm);
		  
void constitutive_update(double * PROPS, MatD2 epsE_t, MatD2 F_t, MatD2 F_tau,
						 MatD2 epsE_tau, MatD2 Sig, MatD4 SpTanMod);
						 
void VMconst(double * PROPS, MatD2 epsE_TR, 
			 MatD2 Dalg, MatD2 Sig, MatD2 epsE_tau);

//void formDsig(MatD2 BeTr, MatD2 kirSig, MatD4 D, MatD2 F_tau, MatD4 SpTanMod, MatD2 Sig);
	  
	  
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
	int nlSdv = *NSVARS;  // number of state variables
	
	MatD2 Fbar_t(3,3);   double detFbar_t;     // previous converged Fbar
	MatD2 Fbar_tau(3,3); double detFbar_tau; // current Fbar
	
	MatD2 F_t(3,3);      double detF_t;       // previous converged def gradient
	MatD2 epsE_t(3,3);                     // previous converged elastic def gradient
	
	MatD2 F_tau(3,3);    double detF_tau;   // current def gradient
	MatD2 epsE_tau(3,3);                   // current elastic def gradient
	
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
			epsE_t = 0.0;
	
        } else {
			
			// ===================================================
			// This is not the first increment; read old values.
            // get previous converged parameters
		    // ===================================================
			jj = 0;
			// elastic deformation gradient
            for ( int i = 0; i < 3; ++i ) {
				for ( int j = 0; j < 3; ++j ) {
					epsE_t(i,j) = svars[jj];
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
		//   epsE_t = previously converged elastic log strain
		//   F_t    = previously converged total def gradient (F-bar previous)
		//   F_tau  = current total def gradient (F-bar current)
		//
		// Output:
		//   epsE_tau = current elastic log strain
		//   T_tau    = current Cauchy stress
		//   SpTanMod = current tangent moduli
		//		
		constitutive_update(PROPS,epsE_t,F_t,F_tau,
							epsE_tau,T_tau,SpTanMod);

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
				svars[jj] = epsE_tau(i,j);
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

void constitutive_update(double * PROPS, MatD2 epsE_t, MatD2 F_t, MatD2 F_tau,
						 MatD2 epsE_tau, MatD2 Sig, MatD4 SpTanMod)
{
	blitz::firstIndex  i;
	blitz::secondIndex j;
	blitz::thirdIndex  k;
	blitz::fourthIndex l;
	
	MatD2 F_t_inv(3,3);  // inverse of total deformation gradient at previous (converged) time, F_t
	MatD2 finc(3,3);     // incremental deformation gradient
	
	VecD  epsEPr_t(3);       // eigenvalues of epsE_t
	MatD2 Nvec_epsE_t(3,3);  // eigenvectors of epsE_t
	
	MatD2 Be_t(3,3);         // left Cauchy-Green tensor at previous converged time, Be_t
	VecD  BePr_t(3);         // eigenvalues of Be_t
	
	MatD2 Be_TR(3,3);        // TRIAL left Cauchy-Green tensor, Be_TR
	VecD  epsBePr_TR(3);     // eigenvalues of Be_TR
	MatD2 Nvec_Be_TR(3,3);   // eigenvectors of Be_TR
	
	MatD2 epsE_TR(3,3);      // TRIAL elastic log strain
	VecD  epsEPr_TR(3);      // eigenvalues of TRIAL elastic log strain
	
	MatD2 kirSig(3,3); 
	MatD2 Dalg(9,9);       // small-strain tangent in MATRIX FORM
	
	// compute incremental deformation gradient
	matlib_inverse(F_t,F_t_inv);
	finc = blitz::sum(F_tau(i,k)*F_t_inv(k,j),k);
	
	// ===========================================================================
	// recover elastic left Cauchy-Green tensor from (half) log elastic strain epsE_t
	// ===========================================================================
	// get eigenvalues and eigenvectors of epsE_t
	solve_eigenvalue_3x3(epsE_t, epsEPr_t, Nvec_epsE_t);

	// compute eigenvalues of Be_t
	for ( int i = 0; i < 3; ++i ) {
		BePr_t(i) = exp(2.0*epsEPr_t(i));
	}
	
	// compute left Cauchy-Green tensor at t 
	Be_t = 0.0;
	for ( int k = 0; k < 3; ++k ) {
		for ( int i = 0; i < 3; ++i ) {
			for ( int j = 0; j < 3; ++j ) {
				Be_t(i,j) += BePr_t(k)*Nvec_epsE_t(k,i)*Nvec_epsE_t(k,j);
			}		
		}
	}
	
	// ===========================================================================
	// compute trial left Cauchy-Green tensor Be_TR = finc * Be_t * finc^T
	// ===========================================================================
	Be_TR = 0.0;
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			for ( int k = 0;  k < 3; ++k ) {
				for ( int l = 0; l < 3; ++l ) {
					Be_TR(i,j) += finc(i,k)*Be_t(k,l)*finc(j,l);
                }
			}
		}
	}
	
	// get eigenvalues and eigenvectors of Be_TR
	solve_eigenvalue_3x3(Be_TR, epsBePr_TR, Nvec_Be_TR);

	// ===========================================================================
	// compute trial elastic (half) log strain tensor
	// ===========================================================================
	// compute principal (half) log strains
	for ( int i = 0; i < 3; ++i ) {
		epsEPr_TR(i) = 0.5*log(epsBePr_TR(i));
	}
	
    // compute epsE_TR 
	epsE_TR = 0.0;
	for ( int k = 0; k < 3; ++k ) {
		for ( int i = 0; i < 3; ++i ) {
			for ( int j = 0; j < 3; ++j ) {
				epsE_TR(i,j) += epsEPr_TR(k)*Nvec_Be_TR(k,i)*Nvec_Be_TR(k,j);
			}		
		}
	}
	
	
	// small-strain return map to get small-strain consistent
    // tangent, Kirchhoff stress, and elastic small-strain tensor
    VMconst(PROPS,epsE_TR,
			Dalg,kirSig,epsE_tau);
                
    // convert to spatial consistent tangent and Cauchy stress
	// note: F_tau is the trial deformation gradient i.e. F_tau = finc*F_t
    //formDsig(Be_TR,kirSig,Dalg,F_tau,SpTanMod,Sig);
    
	
	
}

void VMconst(double * PROPS, MatD2 epsE_TR, 
			 MatD2 Dalg, MatD2 Sig, MatD2 epsE_tau)
{
	// BELOW IS UNCHECKED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// references
	// 1) Ziekienwicz,Taylor FEM 2nd volume, 7th edition - stuff on 9x9 D and size-9 vectors
	// 2) Comp plas, Souza, et al
	//
	double tol = 1.e-6;
	double E,nu,sigmay;
	double Kmod,Gmod;
	double sdevnorm,trSig,j2,f;
	EgVecD bm1; bm1.resize(9); bm1.setZero();
	
	EgVecD epsE9;   epsE9.resize(9); 
	EgVecD epsEtr9; epsEtr9.resize(9); 
	EgVecD Sig9;    Sig9.resize(9);
	EgVecD sdev;    sdev.resize(9); 
	
	
	//EgMatD2 De6;   De6.resize(6,6);   De6.setZero();
	//EgMatD2 Dalg6; Dalg6.resize(6,6); Dalg6.setZero();
	
	EgMatD2 Dalg9; Dalg9.resize(9,9); Dalg9.setZero();
	EgMatD2 De9;   De9.resize(9,9);   De9.setZero();
	EgMatD2 Ce9;   Ce9.resize(9,9);   Ce9.setZero();
	EgMatD2 Id9;   Id9.resize(9,9);   Id9.setZero();
	EgMatD2 Idev;  Idev.resize(9,9);  Idev.setZero();
	EgMatD2 PmatT; PmatT.resize(6,9); PmatT.setZero(); // matrix to convert between 9- and 6-element stuff
	
	// identity tensor size-9 vector form
	bm1(0) = 1.0;
	bm1(1) = 1.0;
	bm1(2) = 1.0;
	
	for ( int i = 0; i < 9; ++i ) {
		Id9(i,i) = 1.0;
	}
	Idev = Id9 - bm1*bm1.transpose()/3.0;
	
	PmatT << 1.0,  0,  0,  0,  0,  0,  0,  0,  0,
		 	   0,1.0,  0,  0,  0,  0,  0,  0,  0,
			   0,  0,1.0,  0,  0,  0,  0,  0,  0,
			   0,  0,  0,0.5,0.5,  0,  0,  0,  0,
			   0,  0,  0,  0,  0,0.5,0.5,  0,  0,
			   0,  0,  0,  0,  0,  0,  0,0.5,0.5;
	
	// Obtain material properties
    E      = PROPS[0];
	nu     = PROPS[1];
	sigmay = PROPS[2];

	Gmod = E/(2.0*(1.0+nu));     // shear modulus
	Kmod = E/(3.0*(1.0-2.0*nu)); // bulk modulus
	
	// convert elastic trial strain tensor to a size-9 vector		 
	epsEtr9(0) = epsE_TR(0,0); // xx  11
	epsEtr9(1) = epsE_TR(1,1); // yy  22
	epsEtr9(2) = epsE_TR(2,2); // zz  33
	epsEtr9(3) = epsE_TR(0,1); // xy  12
	epsEtr9(4) = epsE_TR(1,0); // yx  21
	epsEtr9(5) = epsE_TR(1,2); // yz  23
	epsEtr9(6) = epsE_TR(2,1); // zy  32
	epsEtr9(7) = epsE_TR(2,0); // zx  31
	epsEtr9(8) = epsE_TR(0,2); // zx  13
	
	// initialize elastic strain to trial elastic strain
	for ( int i = 0; i <= 8; ++i ) {
		epsE9(i) = epsEtr9(i);
	}
	
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
	sdev = Sig9-trSig*bm1;
	
	// Compute J2=1/2*s_ij*s_ij
	sdevnorm = sdev.norm();
	j2 = 0.5*sdevnorm*sdevnorm;
	
	// yield function value
	f = sqrt(3.0*j2)/sigmay-1.0;

	if ( f > tol )
	{
		// ===============================
		// PLASTIC UPDATE
		// ===============================
		
		EgVecD resid; resid.resize(10); 
		EgVecD    df; df.resize(9);
		EgMatD2  ddf; ddf.resize(9,9);
		EgMatD2  dgam_ddf_De9; dgam_ddf_De9.resize(9,9);
		
		EgMatD2    A; A.resize(10,10);
		EgVecD    dx; dx.resize(10);
		
		int    itnum, maxit;
		double dgam;
		maxit = 25;
		resid.setZero(); // initialize residual vector
		resid(9)  = f;
		itnum     = 0;     // iteration counter
		dgam      = 0.0;   // plastic multiplier
		
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
		df  = sqrt(3.0/2.0)/sigmay*sdev/sdevnorm;
		ddf = sqrt(3)/(2.0*sigmay*sqrt(j2))*( Idev - sdev*sdev.transpose()/(2.0*j2) ); 
   
		// newton iteration
		while ( ( itnum < maxit ) && ( resid.head(9).norm() > tol ) || ( abs(resid(9))>tol ) )
		{	
			itnum++;
			
			// form Jacobian
			dgam_ddf_De9.setZero();
			for ( int i = 0; i <= 8; ++i ) {
				for ( int j = 0; j <= 8; ++j ) {
					for ( int k = 0; k <= 8; ++k ) {
						dgam_ddf_De9(i,j) += dgam*ddf(i,k)*De9(k,j);
					}
				}
			}
			
			A.setZero();
			for ( int i = 0; i <= 8; ++i ) {
				for ( int j = 0; j <= 8; ++j ) {
					A(i,j) = Id9(i,j) + dgam_ddf_De9(i,j);
				}
			}
			
			for ( int i = 0; i <= 8; ++i ) {
				for ( int j = 0; j <= 8; ++j ) {
					A(9,i) += df(j)*De9(j,i);
				}
				A(i,9) = df(i);	
			}
			
			
			dx = -A.fullPivLu().solve(resid);
			
			epsE9 += dx.head(9); // elems 1 through 9
			dgam  += dx(9);      // elem 10
			
			Sig9 = De9*epsE9;
			
			// deviatoric stress vector
			trSig = (Sig9(0)+Sig9(1)+Sig9(2))/3.0;
			sdev = Sig9-trSig*bm1;
	
			// Compute J2=1/2*s_ij*s_ij
			sdevnorm = sdev.norm();
			j2 = 0.5*sdevnorm*sdevnorm;
	
			// recompute terms in jacobian matrix for next iteration
			df  = sqrt(3.0/2.0)/sigmay*sdev/sdevnorm;
			ddf = sqrt(3.0)/(2.0*sigmay*sqrt(j2))*( Idev - sdev*sdev.transpose()/(2.0*j2) ); 
   
			// recompute residual vector
			for ( int i = 0; i <= 8; ++i ) {
				resid(i) = epsE9(i)-epsEtr9(i) + dgam*df(i);
			}
			resid(9) = sqrt(3.0*j2)/sigmay-1.0;
			
		} // end return mapping 
		
		// consistent tangent operator
		A.setZero();
		for ( int i = 0; i <= 8; ++i ) {
			for ( int j = 0; j <=8 ; ++j ) {
				A(i,j) = Ce9(i,j) + dgam*ddf(i,j);
			}
			A(9,i) = df(i);
			A(i,9) = df(i);			
		}
		
		EgMatD2 B = A.inverse();
		Dalg9 = B.topLeftCorner(9,9); 
		//Dalg6 = PmatT*Dalg9*PmatT.transpose();
		
		for ( int i = 0; i < 9; ++i ) {
			for ( int j = 0; j < 9; ++j ) {
				Dalg(i,j) = Dalg9(i,j);
			}
		}
		
	} else {
		
		// ===============================
		// ELASTIC UPDATE
		// ===============================
		//De6 = PmatT*De9*PmatT.transpose();
		
		for ( int i = 0; i < 9; ++i ) {
			for ( int j = 0; j < 9; ++j ) {
				Dalg(i,j) = De9(i,j);
			}
		}
		

	} // end check yield	
	
	// convert strains from 9 to output form
	epsE_tau(0,0) = epsE9(0); // xx
	epsE_tau(1,1) = epsE9(1); // yy
	epsE_tau(2,2) = epsE9(2); // zz
	epsE_tau(0,1) = epsE9(3); // xy
	epsE_tau(1,0) = epsE9(4); // yx
	epsE_tau(1,2) = epsE9(5); // yz
	epsE_tau(2,1) = epsE9(6); // zy
	epsE_tau(2,0) = epsE9(7); // zx
	epsE_tau(0,2) = epsE9(8); // xz
	
	// convert stresses from 9 to output form
	Sig(0,0) = Sig9(0);
	Sig(1,1) = Sig9(1);
	Sig(2,2) = Sig9(2);
	Sig(0,1) = Sig9(3);
	Sig(1,0) = Sig9(4);
	Sig(1,2) = Sig9(5);
	Sig(2,1) = Sig9(6);
	Sig(2,0) = Sig9(7);
	Sig(0,2) = Sig9(8);
	

}

//void formDsig(MatD2 BeTr, MatD2 kirSig, MatD4 D, MatD2 F_tau, MatD4 SpTanMod, MatD2 Sig)
//{	
//}

void solve_eigenvalue_3x3(MatD2 AMat, VecD BePr, MatD2 Nvec)
{
	
	Eigen::EigenSolver<Eigen::Matrix3d> es;
	Eigen::Matrix3d A;
	
	// assign AMat to A
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			A(i,j) = AMat(i,j);
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
		BePr(i) = EigArray[i].value;
		for ( int j = 0; j < 3; ++j ) {
			Nvec(i,j) = EigArray[i].vector[j];
		}
	}
	
}


// Determinant of a 3x3 matrix
double matlib_determinant(MatD2 A)
{
	double det;

	// 0 1 2 = [0][0] [0][1] [0][2]
	// 3 4 5 = [1][0] [1][1] [1][2]
	// 6 7 8 = [2][0] [2][1] [2][2]
	det = A(0,0)*(A(1,1)*A(2,2)-A(1,2)*A(2,1))
         -A(0,1)*(A(1,0)*A(2,2)-A(1,2)*A(2,0))
         +A(0,2)*(A(1,0)*A(2,1)-A(1,1)*A(2,0));
   
	return det;
}

double matlib_inverse(MatD2 A, MatD2 Ainv)
{
	double det;

	det = matlib_determinant(A);
	
	if (abs(det) < 1.e-20) return 0.e0;

	// 0 1 2 = [0][0] [0][1] [0][2]
	// 3 4 5 = [1][0] [1][1] [1][2]
	// 6 7 8 = [2][0] [2][1] [2][2]

	Ainv(0,0) = ( A(1,1)*A(2,2)-A(1,2)*A(2,1))/det;
	Ainv(0,1) = (-A(0,1)*A(2,2)+A(0,2)*A(2,1))/det;
	Ainv(0,2) = ( A(0,1)*A(1,2)-A(0,2)*A(1,1))/det;
	Ainv(1,0) = (-A(1,0)*A(2,2)+A(1,2)*A(2,0))/det;
	Ainv(1,1) = ( A(0,0)*A(2,2)-A(0,2)*A(2,0))/det;
	Ainv(1,2) = (-A(0,0)*A(1,2)+A(0,2)*A(1,0))/det;
	Ainv(2,0) = ( A(1,0)*A(2,1)-A(1,1)*A(2,0))/det;
	Ainv(2,1) = (-A(0,0)*A(2,1)+A(0,1)*A(2,0))/det;
	Ainv(2,2) = ( A(0,0)*A(1,1)-A(0,1)*A(1,0))/det;

	return det;
}

void onem(MatD2 A)
{
	A = 0.0;
	for ( int i = 0; i < 3; i++ ) A(i,i) = 1.0;
}
