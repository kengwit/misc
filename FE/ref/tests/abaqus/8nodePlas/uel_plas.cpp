/*!************************************************************************
! User element (UEL) for linear elasticity -- plane-strain
!************************************************************************
! Element details:
!************************************************************************
!
! Solution variables (or nodal variables) are the displacements (DOFs 1-2)
!
! Material behavior is linear elasticity 
! 
! This subroutine is for a two-dimensional 8-node isoparametric
!  quadrilateral element as shown below with 4pt (reduced), 9pt (full) 
!
!
!              A eta (=xi_2)
!  8-node      |
!   quad       | Face 3
!        4-----7-----3
!        |     |     |
!        |     |     |
! Face 4 8     ------6---> xi (=xi_1)
!        |           |
!        |           |  Face 2
!        1-----5-----2
!            Face 1
!
!
!************************************************************************
! Usage:
!************************************************************************
!
! User element statement in the input file:
!  *User Element,Nodes=4,Type=U1,Iproperties=1,Properties=2,Coordinates=2,Variables=1,Unsymm
!  1,2
!
! Note: No local state variables are used in this element, thus we may set the above 
!  parameter 'Variables' to any non-zero integer.
!
! In the subroutine UEL, set 'nInt' = number of integration points
!  Options are nInt=4 (reduced), 9 (full integration) 
!
! In the input file, set the parameter pe=1 for plane-strain
!
!     Material Properties Vector
!     --------------------------------------------------------------
!     E  = props(1)  ! Young's modulus
!     nu = props(2)  ! Poisson's ratio
!     pe = jprops(1) ! Plane strain=1, axisymmetric=0
!
!*************************************************************************/


//#include <aba_for_c.h> since 6.13?
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cassert>
using namespace std;
#include "boost/multi_array.hpp"
//
//#include <Eigen/Dense>
//#include <Eigen/Eigenvalues>

/*  
The C++ umat must be compiled externally using a compiler supported by ABAQUS, producing an object file umat_vec.obj. The compiler command is

cl /EHsc /c /Fo uel_plas.cpp /I"C:\cygwin\home\KWL\libraries\eigen-3.2.1" /I"C:\Users\KWL\Desktop\Research\libraries\boost_1_56_0"
cl /EHsc /c /Fo uel_plas.cpp 

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

// column-major for fortran array access
// h = row-height
#define FORT_ACCESS(A,r,c,h) ((A)[(r) + (c)*(h)])

#define NQUAD 4
#define NDM  2
#define NEN  8	
#define NDM1 (NDM+1)
#define NSHP NDM1
#define NSDV 18

void initializeD(double * A, int n, double val);
void compute_shape(int nIntPt, double * shpl, double * jacl, double * xl);
void compute_weight(int nIntPt, double * wl);

double matlib_determinant(double *A);
double matlib_inverse(double *A,double *Ainv);

typedef boost::multi_array<double, 4> Tensor4;
typedef Tensor4::index iT4;

typedef boost::multi_array<double, 2> Tensor2;
typedef Tensor2::index iT2;

typedef boost::multi_array<double, 2> Matrix;
typedef Matrix::index iM2;

void onem(Tensor2 & A);
void Elas(double * props, int nprops, Tensor2 & F, Tensor2 & sig, Tensor4 & SpUUMod);

// note function umat with underscore "uel_" DOES NOT WORK ON MY COMPUTER!
extern "C" void uel(double * RHS, double * AMATRX, double * SVARS, double * ENERGY, int * NDOFEL, int * NRHS, int * NSVARS,
     double * PROPS, int * NPROPS, double * COORDS, int * MCRD, int * NNODE, double * U, double * DU, double * V, double * A, int * JTYPE, double * TIME, double * DTIME,
     int * KSTEP, int * KINC, int * JELEM, double * PARAMS, int * NDLOAD, int * JDLTYP, double * ADLMAG, double * PREDEF,
     int * NPREDF, int * LFLAGS,int * MLVARX,double * DDLMAG,int * MDLOAD,double * PNEWDT,int * JPROPS,int * NJPROP,
     double * PERIOD)
{
	
	assert(JPROPS[0] == 1); // Plane strain=1, axisymmetric=0
	
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
    //  deformation problem
    if(LFLAGS[1]==1) {	
        // LFLAGS[1]=0 -> small disp.
        // LFLAGS[1]=1 -> large disp.
        cout << "Abaqus thinks you are doing" << endl;
        cout << "a large displacement analysis" << endl;
        cout << "go in and remove nlgeom=yes" << LFLAGS[0] << endl;
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
	
	double E,nu,Energy;
	int nrhs,ndofel; 
	int counter;
	int skip1,skip2,skip3;			
	int NDM2 = NDM*NDM;
	
	
	vector<double> shp0;
	vector<double> shpC0; // shape functions at centroid, reference & current
	vector<double> detMapJ0;
	vector<double> detMapJ0C; // determinant jacobians at centroid, reference & current
	vector<double> shp;
	vector<double> shpC;   // shape functions at normal gauss points, reference & current
	vector<double> detMapJ;
	vector<double> detMapJC; // determinant jacobians at normal gauss points, reference & current
	vector<double> w;
	vector<double> u;
	vector<double> du;
	vector<double> coords;
	vector<double> coordsC; // coords = reference, coordsC = current
	vector<double> Ru,Kuu;
	
	// general 3D
	vector<double> body; 
	Tensor2   Iden(boost::extents[3][3]);
	Tensor2  F_old(boost::extents[3][3]);
	Tensor2 Fe_old(boost::extents[3][3]);
	Tensor2 Fe_tau(boost::extents[3][3]);
	Tensor2 Fc_tau(boost::extents[3][3]); double detFc_tau;
	Tensor2  F_tau(boost::extents[3][3]); double detF_tau;
	Tensor2  T_tau(boost::extents[3][3]);
	Tensor4 SpUUMod(boost::extents[3][3][3][3]); 
	
	// specific to plane-strain
	vector<double> BodyForceRes,Smat;
	Matrix Bmat(boost::extents[3][NDM*NEN]);
	Matrix Gmat(boost::extents[4][NDM*NEN]);
	Matrix G0mat(boost::extents[4][NDM*NEN]);
	Matrix Amat(boost::extents[NDM2][NDM2]);
	Matrix Qmat(boost::extents[NDM2][NDM2]);
	cout << "here1\n";
	
	
	ndofel = *NDOFEL; // this is TOTAL dofs (2*8 = 16; 16 displacement dofs) per element
	nrhs   = *NRHS;
	
	// only 1 integration point (at centroid)
	shp0.resize(NSHP*NEN); 
	shpC0.resize(NSHP*NEN); 
	detMapJ0.resize(1); 
	detMapJ0C.resize(1);      
	
	// NQUAD integration points
	shp.resize(NQUAD*NSHP*NEN); 
	shpC.resize(NQUAD*NSHP*NEN); 
	detMapJ.resize(NQUAD); 
	detMapJC.resize(NQUAD);           
	w.resize(NQUAD);
	
	u.resize(NDM*NEN);
	du.resize(NDM*NEN);
	coords.resize(NDM*NEN);
	coordsC.resize(NDM*NEN);
	Ru.resize(NDM*NEN);
	Kuu.resize(NDM*NEN*NDM*NEN);
	
	cout << "here2\n";
	
	// general 3D
	body.resize(3);
	
	// specific to plane-strain
	BodyForceRes.resize(NDM*NEN);
	Smat.resize(3);
	 
	// Identity tensor
    onem(Iden);
	
	// Initialize the residual and tangent matrices to zero.
    initializeD(&Ru[0],NDM*NEN,0.0);
    initializeD(&Kuu[0],NDM*NEN*NDM*NEN,0.0);
    *ENERGY=0.0;

	cout << "here4\n";
	
	// NOTE: COORDS(MCRD,NNODE)
	// COORDS = [0 2 4 6
	//           1 3 5 7 ...]
	// coords[0] = COORDS(1,1)
	// coords[1] = COORDS(2,1)
	// coords[2] = COORDS(1,2)
	// coords[3] = COORDS(2,2)
	// coords[4] = COORDS(1,3)
	// coords[5] = COORDS(2,3)
	
	counter = 0;
	for ( int i = 0; i < NEN; ++i ) {
		for ( int j = 0; j < NDM; ++j ) {
			u[i*NDM+j]  = U[counter];
            du[i*NDM+j] = DU[counter];
			counter++;
			
		}
	}
	
	
	for ( int i = 0; i < NEN; ++i ) {
		for ( int j = 0; j < NDM; ++j ) {
			coords[i*NDM+j] = FORT_ACCESS(COORDS,j,i,NDM);
			coordsC[i*NDM+j] = coords[i*NDM+j];// + u[i*NDM+j]			
		}
	}
	

    for ( int i = 0; i < NEN; ++i ) {
		for ( int j = 0; j < NDM; ++j ) {
			cout << "coords[" << i*NDM+j << "] = " << coords[i*NDM+j] << "  ";
			
		}
		cout << endl;
	}
	
    for ( int i = 0; i < NEN; ++i ) {
		for ( int j = 0; j < NDM; ++j ) {
			cout << "u[" << i*NDM+j << "] = " << u[i*NDM+j] << "  ";
			
		}
		cout << endl;
	}
		
	for ( int i = 0; i < NEN; ++i ) {
		for ( int j = 0; j < NDM; ++j ) {
			cout << "du[" << i*NDM+j << "] = " << du[i*NDM+j] << "  ";
			
		}
		cout << endl;
	}
	
	
	// ----------------------------------------------------------------
	// compute shape functions
	if(NEN == 8) {
		compute_shape(    1, &shp0[0],  &detMapJ0[0],  &coords[0]); // at centroid, reference, for F-bar method
		compute_shape(    1,&shpC0[0], &detMapJ0C[0], &coordsC[0]); // at centroid, current, for F-bar method
        compute_shape(NQUAD,  &shp[0],   &detMapJ[0],  &coords[0]); // at all gauss points, reference
        compute_shape(NQUAD, &shpC[0],  &detMapJC[0], &coordsC[0]); // at all gauss points, current
		compute_weight(NQUAD,&w[0]);
		
    } else {
		cout << "Incorrect number of nodes: nNode.ne.8" << endl;
        exit(1);
	}
	// compute Fc_tau (at centroid) for F-bar method
	// only for 4-noded fullly integrated quad
	if ((NEN==4)&&(NQUAD==4)&&(LFLAGS[1]==1)) 
	{
		cout << "F-BAR\n";
        onem(Fc_tau);
    	for ( int i=0; i < NDM; ++i )
			for ( int j=0; j < NDM; ++j )
				for ( int k=0; k < NEN; ++k )
					Fc_tau[i][j] += u[k*NDM+i]*shp0[k*NDM1+j];
		
		detFc_tau = matlib_determinant(Fc_tau.data());
	}
	
	// ================================================
	// Loop over integration points
    // ================================================
	
	int jj_sv; // to access current block of state variables
    double * svars  = SVARS; // current block of state variables  
    double * sh0_gp  = &shp0[0];
	double * shC0_gp = &shpC0[0];
	
	double * sh_gp  = &shp[0];
	double * j_gp   = &detMapJ[0];
	double * shC_gp = &shpC[0];
	double * jC_gp  = &detMapJC[0];
	
	for ( int ip = 0; ip < NQUAD; ++ip )
	{
		// ================================================
		// Obtain state variables from previous increment
        // ================================================
		if((*KINC <= 1)&&(*KSTEP == 1)) {
            
            // This is the first increment of the first step.
            //  Give initial conditions.
			onem(F_old);
			onem(Fe_old);
	
        } else {
			
			// This is not the first increment; read old values.
            // get previous converged parameters
		    jj_sv = 0;
			// total deformation gradient
            for ( int i = 0; i < 3; ++i ) {
				for ( int j = 0; j < 3; ++j ) {
					F_old[i][j] = svars[jj_sv];
					jj_sv++;
				}
			}
			// elastic deformation gradient
            for ( int i = 0; i < 3; ++i ) {
				for ( int j = 0; j < 3; ++j ) {
					Fe_old[i][j] = svars[jj_sv];
					jj_sv++;
				}
			}
			
        }
	
		// ==================================================================
		// Obtain the deformation gradient at this integration point.
        //  The subscript tau denotes the time at the end of the increment.
        //
        onem(F_tau);
        for ( int i=0; i < NDM; ++i )
			for ( int j=0; j < NDM; ++j )
				for ( int k=0; k < NEN; ++k )
					F_tau[i][j] += u[k*NDM+i]*sh_gp[k*NDM1+j];
				
		// If large deformation:
        // Modify the deformation gradient for the `F-bar' method
        //  only when using the 8 node fully integrated linear
        //  element, do not use the `F-bar' method for reduced element
        // only for 4-noded fullly integrated quad
		if ((NEN==4)&&(NQUAD==4)&&(LFLAGS[1]==1)) 
		{
			cout << "Fbar not yet implemented\n";
            //detF_tau = matlib_determinant(F_tau);
			//double fac = pow((detFc_tau/detF_tau),(1.0/3.0));
			//for ( int i = 0; i < 9; ++i ) F_tau[i] = fac*F_tau[i];
			
        }
		
		cout << "F_tau =============\n";
		for ( int i=0; i < 3; ++i ) {
			for ( int j=0; j < 3; ++j ) {
				cout << F_tau[i][j] << " ";
			}
			cout << endl;
		}
		
		//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        //
        // Perform the constitutive update at this integ. point
        // - returns spatial tangent moduli
        Elas(PROPS,*NPROPS,F_tau,T_tau,SpUUMod);
        //
        //
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


        // Compute/update the displacement residual vector
        if (JPROPS[0] == 1) 
		{
			// this is plane strain
            Smat[0] = T_tau[0][0];
            Smat[1] = T_tau[1][1];
            Smat[2] = T_tau[0][1];
			
			initializeD(Bmat.data(),3*NDM*NEN,0.0);
            for ( int k=0; k < NEN; ++k )
			{
				skip1 = NDM*k;
				skip2 = k*NDM1;
				Bmat[0][skip1  ] = shC_gp[skip2  ]; // dNk/dx
				Bmat[1][skip1+1] = shC_gp[skip2+1]; // dNk/dy
				Bmat[2][skip1  ] = shC_gp[skip2+1]; // dNk/dy
				Bmat[2][skip1+1] = shC_gp[skip2  ]; // dNk/dx
            }

            //body = zero ! The body force vector may be specified here
            
            //BodyForceRes = zero
            //do k=1,nNode
            //   BodyForceRes(1+nDim*(k-1),1) = sh(k)*body(1)
            //   BodyForceRes(2+nDim*(k-1),1) = sh(k)*body(2)
            //enddo
            
			// Bmat
			// = [dNk/dx      0]
			//   [     0 dNk/dy]
			//   [dNk/dy dNk/dx]
			
			// Bmat^T
			// = [dNk/dx      0 dNk/dy]
			//   [     0 dNk/dy dNk/dx]
			for ( int k=0; k < NEN*NDM; ++k )
			{
				for ( int j=0; j < 3; ++j )
				{
					Ru[k] = Ru[k] - Bmat[j][k]*Smat[j]*jC_gp[ip]*w[ip];
					//+ BodyForceRes*detmapJC*w(intpt)*AR
				}
			}
            
        } else {
            // big problem
            cout << "not plane strain" << endl;
			exit(1);
            
        }


        // Compute/update the displacement tangent matrix
        //
        if (JPROPS[0] == 1) // this is plane strain
		{	
            
            // Note: Gmat is 4 by nDim*nNode
            initializeD(Gmat.data(),4*NDM*NEN,0.0);
            for ( int k=0; k < NEN; ++k )
			{
				skip1 = NDM*k;
				skip2 = k*NDM1;
				Gmat[0][skip1  ] = shC_gp[skip2  ]; // dNk/dx
				Gmat[1][skip1+1] = shC_gp[skip2  ]; // dNk/dx
				Gmat[2][skip1  ] = shC_gp[skip2+1]; // dNk/dy
				Gmat[3][skip1+1] = shC_gp[skip2+1]; // dNk/dy	
            }

            // Note: G0mat is 4 by nDim*nNode
            initializeD(G0mat.data(),4*NDM*NEN,0.0);
            for ( int k=0; k < NEN; ++k )
			{
				skip1 = NDM*k;
				skip2 = k*NDM1;
				G0mat[0][skip1  ] = shC0_gp[skip2  ]; // dNk/dx
				G0mat[1][skip1+1] = shC0_gp[skip2  ]; // dNk/dx
				G0mat[2][skip1  ] = shC0_gp[skip2+1]; // dNk/dy
				G0mat[3][skip1+1] = shC0_gp[skip2+1]; // dNk/dy	
			}
			
            initializeD(Amat.data(),NDM2*NDM2,0.0);
			
			// ((i*n2+j)*n3+k)*n4+l
			// = i*n2*n3*n4 + j*n3*n4 + k*n4 + l
			// = 27*i + 9*j + 3*k + l
			Amat[0][0] = SpUUMod[0][0][0][0]; //(1,1,1,1)
            Amat[0][1] = SpUUMod[0][0][1][0]; //(1,1,2,1)
            Amat[0][2] = SpUUMod[0][0][0][1]; //(1,1,1,2)
            Amat[0][3] = SpUUMod[0][0][1][1]; //(1,1,2,2) 
            Amat[1][0] = SpUUMod[1][0][0][0]; //(2,1,1,1) 
            Amat[1][1] = SpUUMod[1][0][1][0]; //(2,1,2,1) 
            Amat[1][2] = SpUUMod[1][0][0][1]; //(2,1,1,2) 
            Amat[1][3] = SpUUMod[1][0][1][1]; //(2,1,2,2) 
            Amat[2][0] = SpUUMod[0][1][0][0]; //(1,2,1,1)
            Amat[2][1] = SpUUMod[0][1][1][0]; //(1,2,2,1)
            Amat[2][2] = SpUUMod[0][1][0][1]; //(1,2,1,2)
            Amat[2][3] = SpUUMod[0][1][1][1]; //(1,2,2,2)
            Amat[3][0] = SpUUMod[1][1][0][0]; //(2,2,1,1)
            Amat[3][1] = SpUUMod[1][1][1][0]; //(2,2,2,1)
            Amat[3][2] = SpUUMod[1][1][0][1]; //(2,2,1,2)
            Amat[3][3] = SpUUMod[1][1][1][1]; //(2,2,2,2)
            /*!
            Qmat = zero
            Qmat(1,1) = half*(Amat(1,1)+Amat(1,4)) - half*T_tau(1,1)
            Qmat(2,1) = half*(Amat(2,1)+Amat(2,4)) - half*T_tau(1,2)
            Qmat(3,1) = half*(Amat(3,1)+Amat(3,4)) - half*T_tau(1,2)
            Qmat(4,1) = half*(Amat(4,1)+Amat(4,4)) - half*T_tau(2,2)
            Qmat(1,4) = half*(Amat(1,1)+Amat(1,4)) - half*T_tau(1,1)
            Qmat(2,4) = half*(Amat(2,1)+Amat(2,4)) - half*T_tau(1,2)
            Qmat(3,4) = half*(Amat(3,1)+Amat(3,4)) - half*T_tau(1,2)
            Qmat(4,4) = half*(Amat(4,1)+Amat(4,4)) - half*T_tau(2,2)
            !
            if((nNode.eq.4).and.(nInt.eq.4)) then
               !
               ! This is the tangent using the F-bar method with the
               !  4 node fully integrated linear element
               !
               Kuu = Kuu
     +              + matmul(matmul(transpose(Gmat),Amat),Gmat)
     +              *detMapJC*w(intpt)*AR
!     +              +matmul(transpose(Gmat),matmul(Qmat,(G0mat - Gmat)))
!     +              *detMapJC*w(intpt)*AR
            else
               !
               ! This is the tangent not using the F-bar method with all
               !  other elements
               !
			   ! big problem
               write(*,*) 'F-bar tangent term not implemented'
               call xit 
			   
            endif
            !*/
        } else {
			// big problem
            cout << "not plane strain" << endl;
			exit(1);
		 }
	
		// ================================================
		// Store current state (Abaqus will 
		// automatically perform update at convergence)
		// ================================================
		jj_sv = 0;
		// total deformation gradient
        for ( int i = 0; i < 3; ++i ) {
			for ( int j = 0; j < 3; ++j ) {
				svars[jj_sv] = F_tau[i][j];
				jj_sv++;
			}
		}
		
		// elastic deformation gradient
        for ( int i = 0; i < 3; ++i ) {
			for ( int j = 0; j < 3; ++j ) {
				svars[jj_sv] = Fe_tau[i][j];
				jj_sv++;
			}
		}
         
		// ================================================
		// Take strides
        // ================================================
		sh_gp += NSHP*NEN;
		shC_gp += NSHP*NEN;
		
		svars += NSDV; 
		
	}
	

}

void compute_weight(int nIntPt, double * wl)
{
	assert ( nIntPt == 1 || nIntPt == 4 || nIntPt == 9 );
	
	
	if ( nIntPt == 1 ) {
		
		wl[0] = 2*2;
		
	} else if ( nIntPt == 4 ) {
	
		wl[0] = 1.0*1.0;
        wl[1] = 1.0*1.0;
		wl[2] = 1.0*1.0; 
	    wl[3] = 1.0*1.0;
		
	} else if ( nIntPt == 9 ) {
		
		double w1 = 5./9.; // corresponding to pt
		double w2 = 8./9.; // corresponding to 0.0
	
		wl[0] = w1*w1;
        wl[1] = w2*w1;
		wl[2] = w1*w1; 
	    wl[3] = w1*w2;
		wl[4] = w2*w2;
		wl[5] = w1*w2;
		wl[6] = w1*w1;
		wl[7] = w2*w1;
		wl[8] = w1*w1;
	
	} else {
		cout << "compute_weight: shouldn't be here\n";	
	}
	
	
}


void compute_shape(int nIntPt, double * shpl, double * jacl, double * xl)
{
	assert ( nIntPt == 1 || nIntPt == 4 || nIntPt == 9 );
	
	double pt; 
	
	
	double * s;
	
	s = new double[nIntPt*NDM];
	
	if ( nIntPt == 1 )
	{
		s[0] = 0.0;
		s[1] = 0.0;
		
	} else if ( nIntPt == 4 ) {
		
		// location of gauss points
		// 3 4
		// 1 2
		
		pt = 1./sqrt(3.);
		s[0] = -pt;
		s[1] = -pt;
		
		s[2] = pt;
		s[3] = -pt;
		
		s[4] = pt;
		s[5] = pt;
		
		s[6] = -pt;
		s[7] = pt;
		
	} else if ( nIntPt == 9 ) {
		
		// location of gauss points
		// 7 8 9
		// 4 5 6
		// 1 2 3
	
		pt = sqrt(0.6);
		s[ 0] = -pt; 
		s[ 1] = -pt;
		s[ 2] = 0.;
		s[ 3] = -pt;
		s[ 4] = pt;
		s[ 5] = -pt;
		s[ 6] = -pt;
		s[ 7] = 0.;
		s[ 8] = 0.;
		s[ 9] = 0.;
		s[10] = pt;
		s[11] = 0.;
		s[12] = -pt;
		s[13] = pt; 
	    s[14] = 0.;
		s[15] = pt;
		s[16] = pt;
		s[17] = pt;

	} else {
		cout << "compute_shape: shouldn't be here\n";
	}
		
	
	vector<double> dnds(NDM*NEN); // temp
	vector<double> xjac(NDM*NDM); // temp
	vector<double> dxds(NDM*NDM); // temp
	double xsj;           		  // temp
	
	int shapeinc = NSHP*NEN;   // (ndm+1)*nen = 3*8 = 24;
	
	double * shp_gp = shpl;  // shape block at current gauss pt
	double * jac_gp = jacl;  // jacobian at current gauss pt
	
	double * s_gp   = s;     // current gauss pt coordinate
		
	//double temp=0;
	for (int lquad = 0; lquad < nIntPt; lquad++)
	{
		// shape functions itself 
		shp_gp[2]  = 0.25*(1.0-s_gp[0])*(1.0-s_gp[1])*(-s_gp[0]-s_gp[1]-1.0);
		shp_gp[5]  = 0.25*(1.0+s_gp[0])*(1.0-s_gp[1])*( s_gp[0]-s_gp[1]-1.0);
		shp_gp[8]  = 0.25*(1.0+s_gp[0])*(1.0+s_gp[1])*( s_gp[0]+s_gp[1]-1.0);
		shp_gp[11] = 0.25*(1.0-s_gp[0])*(1.0+s_gp[1])*(-s_gp[0]+s_gp[1]-1.0);
		
		shp_gp[14] = 0.5*(1.0-s_gp[1])*(1.0+s_gp[0])*(1.0-s_gp[0]);
		shp_gp[17] = 0.5*(1.0-s_gp[1])*(1.0+s_gp[1])*(1.0+s_gp[0]);
		shp_gp[20] = 0.5*(1.0+s_gp[1])*(1.0+s_gp[0])*(1.0-s_gp[0]);
		shp_gp[23] = 0.5*(1.0-s_gp[0])*(1.0+s_gp[1])*(1.0-s_gp[1]);
	
		// ===============================
		// Compute shape function derivatives wrt natural coordinates 
		//  - this is an intermediate stage
		// ===============================
		// shape function derivatives (d/ds)
		// For e.g.:
		// dN/ds = [dN1/dr ... dNnen/dr;
		// 		    dN1/ds ... dNnen/ds]; ndm*nen
		// remember that u = 1-r-s-t
		dnds[ 0] =   s_gp[0]*( 0.5 - 0.5*s_gp[1] ) + 0.25*s_gp[1] - 0.25*s_gp[1]*s_gp[1]; // dN1/ds
		dnds[ 1] =  -0.25*s_gp[0]*s_gp[0] + s_gp[0]*(0.25 - 0.5*s_gp[1]) + 0.5*s_gp[1];   // dN1/dt

		dnds[ 2] =  s_gp[0]*(0.5 - 0.5*s_gp[1]) - 0.25*s_gp[1] + 0.25*s_gp[1]*s_gp[1]; // dN2/ds
		dnds[ 3] =  -0.25*s_gp[0]*s_gp[0] + s_gp[0]*(-0.25 + 0.5*s_gp[1]) + 0.5*s_gp[1]; // dN2/dt

		dnds[ 4] = 0.25*(1. + s_gp[1])*(2.*s_gp[0] + 1.*s_gp[1]); // dN3/ds
		dnds[ 5] = 0.25*(1. + s_gp[0])*(1.*s_gp[0] + 2.*s_gp[1]); // dN3/dt

		dnds[ 6] = s_gp[0]*(0.5 + 0.5*s_gp[1]) - 0.25*s_gp[1] - 0.25*s_gp[1]*s_gp[1];
		dnds[ 7] = 0.25*s_gp[0]*s_gp[0] + s_gp[0]*(-0.25 - 0.5*s_gp[1]) + 0.5*s_gp[1];
		
		dnds[ 8] = s_gp[0]*(-1. + 1.*s_gp[1]);
		dnds[ 9] = -0.5 + 0.5*s_gp[0]*s_gp[0];
		
		dnds[10] = 0.5 - 0.5*s_gp[1]*s_gp[1];
		dnds[11] = (-1. - 1.*s_gp[0])*s_gp[1];
		
		dnds[12] = s_gp[0]*(-1. - 1.*s_gp[1]);
		dnds[13] = 0.5 - 0.5*s_gp[0]*s_gp[0];
		
		dnds[14] = -0.5 + 0.5*s_gp[1]*s_gp[1];
		dnds[15] = (-1. + 1.*s_gp[0])*s_gp[1];
		
		// ==============================
		// Compute Jacobian matrix
		// ==============================
		// 3D dxds = [dx/dr dy/dr dz/dr; = [0 1 2 
		//            dx/ds dy/ds dz/ds;    3 4 5
		//			  dx/dt dy/dt dz/dt]    6 7 8];
		//
		// 2D dxds = [dx/ds dy/ds;  = [0 1
		//            dx/dt dy/dt];    2 3];
		//
		// 1D dxds = [dx/ds] 
		// 
		for (int i=0; i < NDM*NDM; i++) {
			dxds[i] = 0.0;
		}
		
		int skip,skip1;
		for (int i=0; i < NEN; i++) {
			skip=i*NDM;
			
			dxds[0] += dnds[skip  ]*xl[skip  ]; // dx/ds = dN1/ds*x1 + ... + dN8/ds*x3 
			dxds[1] += dnds[skip  ]*xl[skip+1]; // dy/ds = dN1/ds*y1 + ... + dN8/ds*y3 
			
			dxds[2] += dnds[skip+1]*xl[skip  ]; // dx/dt = dN1/dt*x1 + ... + dN8/dt*x3 
			dxds[3] += dnds[skip+1]*xl[skip+1]; // dy/dt = dN1/dt*y1 + ... + dN8/dt*y3 
			
		}

		// ==============================
		// Compute determinant of Jacobian matrix
		// ==============================
		xsj = dxds[0]*dxds[3]-dxds[2]*dxds[1];

		// ==============================
		// Compute INVERSE of Jacobian matrix
		// ==============================
		// note: inverse of Jac is NOT just the Jac with terms inversed!
		// 3D xjac = [dr/dx ds/dx dt/dx; = [0 1 2
		//            dr/dy ds/dy dt/dy     3 4 5
		//            dr/dz ds/dz dt/dz];   6 7 8];
		//
		// 2D xjac = [ds/dx dt/dx;  = [0 1
		//            ds/dy dt/dy];    2 3];
		//
		// 1D xjac = [ds/dx];
		//		
		// below is xjac = cofactor(dxds)^Transpose / det(dxds)
		xjac[0] = dxds[3]/xsj; // ds/dx
		xjac[1] =-dxds[1]/xsj; // dt/dx
		xjac[2] =-dxds[2]/xsj; // ds/dy
		xjac[3] = dxds[0]/xsj; // dt/dy
		
		// ==============================
		// store determinant of Jacobian 
		// matrix times integration weight
		// ==============================
		jac_gp[lquad] = xsj; // xsj=ratio of real space area to reference area (4 units) e.g. quarter for unit square in real space
    	
		// ==============================
		// compute and store shape function 
		// derivatives wrt to real coordinates
		// ==============================
		for ( int i = 0; i < NEN; ++i ) {
			skip  = i*NDM;
			skip1 = i*NDM1;
			
			// just think of Jac^inv * {dNi/ds}
			//                         {dNi/dt} 	
			shp_gp[skip1  ] = xjac[0]*dnds[skip]+xjac[1]*dnds[skip+1]; // dNi/dx = dNi/ds*ds/dx + dNi/dt*dt/dx
			shp_gp[skip1+1] = xjac[2]*dnds[skip]+xjac[3]*dnds[skip+1]; // dNi/dy = dNi/ds*ds/dy + dNi/dt*dt/dy
		}

		
		// ==============================
		//  take stride 
		// ==============================
		shp_gp += shapeinc;
		s_gp   += 2;
		
		 
	}
	
	delete [] s; s = NULL;
	
}


void initializeD(double * A, int n, double val)
{
	for (int i = 0; i < n; ++i ) A[i] = val;
}

void onem(Tensor2 & A)
{
	for (int i = 0; i < 3; ++i )
		for (int j = 0; j < 3; ++j )
			A[i][j] = 0.0;
	
	for (int i = 0; i < 3; ++i ) A[i][i] = 1.0;
	
}



//-----CONVENIENCE FUNCTIONS-------------------------------------------------

// Determinant of a 3x3 matrix
double matlib_determinant(double *A)
{
  double det;

  det = A[0]*(A[4]*A[8]-A[5]*A[7])
       -A[1]*(A[3]*A[8]-A[5]*A[6])
       +A[2]*(A[3]*A[7]-A[4]*A[6]);

  return det;
}

// Inversion of a 3x3 matrix, which also
// returns the determinant as a by-product.
// If the determinant is close to zero,
// no inversion is performed.
double matlib_inverse(double *A,double *Ainv)
{
  double det,detinv;

  det = matlib_determinant(A);
  if (fabs(det) < 1.e-20) return 0.e0;

  detinv = 1./det;
  Ainv[0] = detinv*( A[4]*A[8]-A[5]*A[7]);
  Ainv[1] = detinv*(-A[1]*A[8]+A[2]*A[7]);
  Ainv[2] = detinv*( A[1]*A[5]-A[2]*A[4]);
  Ainv[3] = detinv*(-A[3]*A[8]+A[5]*A[6]);
  Ainv[4] = detinv*( A[0]*A[8]-A[2]*A[6]);
  Ainv[5] = detinv*(-A[0]*A[5]+A[2]*A[3]);
  Ainv[6] = detinv*( A[3]*A[7]-A[4]*A[6]);
  Ainv[7] = detinv*(-A[0]*A[7]+A[1]*A[6]);
  Ainv[8] = detinv*( A[0]*A[4]-A[1]*A[3]);

  return det;
}


void Elas(double * props, int nprops, Tensor2 & F1, Tensor2 & sig, Tensor4 & SpUUMod)
{
	for ( int i = 0; i < 3; ++i )
	for ( int j = 0; j < 3; ++j )
	for ( int k = 0; k < 3; ++k )
	for ( int l = 0; l < 3; ++l )
		SpUUMod[i][j][k][l]=0.;	
}
