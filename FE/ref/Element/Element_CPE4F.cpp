#include "Element_CPE4F.h"

// =================================================================================
//                          2D 4-NODE F-BAR PLANE STRAIN ELEMENT 
// =================================================================================
Element_CPE4F::Element_CPE4F(unsigned int id) : Element(id)
{ 
	// cout << "Destroy CPE4F element " << el_index << endl; 
	element_formulation = FBAR;
	section_type = PLANE_STRAIN;
	nen          = 4;                   // nodes per element
	ndm          = 2;                   // number of dimensions
	ndf          = 2;                   // number of dof per node
	nquad        = 4;                   // number of quadrature points
	nsdm         = ndm*ndf;             // strain dimension = 2*2 = 4
	nst          = ndf*nen;             // number of dof per element
	nvoigt       = ndm*(ndm+1)/2;       // number of stress terms in voigt notation
	
	assert(nvoigt == 3);
	assert(fem::formulation == fem::updated_lagrangian);
	allocate_data();
	
	
} 

Element_CPE4F::~Element_CPE4F() 
{ 
	// cout << "Destroy CPE4F element " << el_index << endl; 
	clear_data();
}
 


void Element_CPE4F::compute_strains( int nintpts, MatrixD3 & F, MatrixD3 & dshp, MatrixD2 & x, MatrixD2 & u, MatrixD2 & du, MatrixD2 & q )
{
	// zero out first
	std::fill( F.origin(), F.origin() + F.num_elements(), 0.0 );
	
	if ( fem::formulation == fem::updated_lagrangian )
	{
		// compute deformation gradients FiJ at time n+1
		for ( int lquad = 0; lquad < nintpts; lquad++ )
		{
				
			for(int i = 0; i < ndf; i++) {
				for(int J = 0; J < ndm; J++) {
					for(int a = 0; a < nen; a++) {
						F[lquad][i][J]+=(x[a][i]+u[a][i])*dshp[lquad][a][J]; 
					}
				}    
			}
			
			F[lquad][2][2] = 1.0; // plane strain	condition
					
		}
	} else {
		cout << "CPE4F - small strains not applicable to this element. Exiting ...\n";
		exit(1);
	}
	
}	
	
void Element_CPE4F::weights(int nintpts, VectorD & w)
{
	assert ( nintpts == 1 || nintpts == 4 );
	
	if ( nintpts == 1 ) {
		
		w[0] = 4.0;//2.0*2.0;
		
	} else if ( nintpts == 4 ) {
	
		w[0] = 1.0;//1.0*1.0;
        w[1] = 1.0;//1.0*1.0;
		w[2] = 1.0;//1.0*1.0; 
	    w[3] = 1.0;//1.0*1.0;
		
	} else {
		cout << "weight CPE4F: shouldn't be here\n";	
	}
	
}

void Element_CPE4F::compute_Bmat(int nintpts, MatrixD3 & B, MatrixD3 & dsh)
{
	std::fill(  B.origin(),  B.origin() +  B.num_elements(), 0.0 );
	
	int skip;
		
	// see Peric complas book Eq. 4.30 pg 89 - we follow convention of this book
	for (int lquad = 0; lquad < nintpts; lquad++)
	{
		for (int a = 0; a < nen; a++)
		{
			skip = a*ndf;
			B[lquad][0][skip  ] = dsh[lquad][a][0]; // d/dx -> 11
			B[lquad][1][skip+1] = dsh[lquad][a][1]; // d/dy -> 22
			B[lquad][2][skip  ] = dsh[lquad][a][1]; // d/dy -> 12
			B[lquad][2][skip+1] = dsh[lquad][a][0]; // d/dx -> 21
		}		
	}

			
}

void Element_CPE4F::compute_Gmat(int nintpts, MatrixD3 & G, MatrixD3 & dsh)
{
	std::fill(  G.origin(),  G.origin() +  G.num_elements(), 0.0 );
	
	int skip;
	
	// see Peric complas book Eq. 4.94 - we follow convention of this book
	for (int lquad = 0; lquad < nintpts; lquad++)
	{
		for (int a = 0; a < nen; a++)
		{
			skip = a*ndf;
			G[lquad][0][skip  ] = dsh[lquad][a][0]; // d/dx -> 11
			G[lquad][1][skip+1] = dsh[lquad][a][0]; // d/dx -> 21
			G[lquad][2][skip  ] = dsh[lquad][a][1]; // d/dy -> 12 
			G[lquad][3][skip+1] = dsh[lquad][a][1]; // d/dy -> 22
		}		
	}
	
		
}

void Element_CPE4F::process_FBarKinematics(ViewD2 & Fc_t, ViewD2 & Fc_tau, ViewD2 & F_t, ViewD2 & F_tau)
{
	// 2D plane-strain implementation detF
	double detFc_t   = Fc_t[0][0]*Fc_t[1][1] - Fc_t[0][1]*Fc_t[1][0];
	double detFc_tau = Fc_tau[0][0]*Fc_tau[1][1] - Fc_tau[0][1]*Fc_tau[1][0];


	// Modify the deformation gradient for the `F-bar' method
	// only when using the 4 node fully integrated linear
	// element, do not use the `F-bar' method for any other element
	//
	//  2D plane-strain implementation
	double detF_t   = F_t[0][0]*F_t[1][1] - F_t[0][1]*F_t[1][0];
	double detF_tau = F_tau[0][0]*F_tau[1][1] - F_tau[0][1]*F_tau[1][0];
	double fac_t   = pow((detFc_t/detF_t),0.5);
	double fac_tau = pow((detFc_tau/detF_tau),0.5);
	
	// note below the factor does not multiply F33 (which remains at 1.0)
	for ( int i = 0; i < ndm; ++i ) {
		for ( int j = 0; j < ndm; ++j ) { 
			F_tau[i][j] = fac_tau*F_tau[i][j];
			F_t[i][j] = fac_t*F_t[i][j];
		}
	}
}

void Element_CPE4F::compute_Qmat(ViewD2 & Qmat, ViewD2 & Amat, ViewD1 & Svec)
{
	//cout << "Qmat.num_elements() = " << Qmat.num_elements() << endl; 
	fill(Qmat.origin(),Qmat.origin()+Qmat.num_elements(),0.0);
	Qmat[0][0] = 0.5*( Amat[0][0]+Amat[0][3] ) - 0.5*Svec[0]; //11
	Qmat[1][0] = 0.5*( Amat[1][0]+Amat[1][3] ) - 0.5*Svec[2]; //12
	Qmat[2][0] = 0.5*( Amat[2][0]+Amat[2][3] ) - 0.5*Svec[2]; //12
	Qmat[3][0] = 0.5*( Amat[3][0]+Amat[3][3] ) - 0.5*Svec[1]; //22
	
	Qmat[0][3] = 0.5*( Amat[0][0]+Amat[0][3] ) - 0.5*Svec[0];
	Qmat[1][3] = 0.5*( Amat[1][0]+Amat[1][3] ) - 0.5*Svec[2];
	Qmat[2][3] = 0.5*( Amat[2][0]+Amat[2][3] ) - 0.5*Svec[2];
	Qmat[3][3] = 0.5*( Amat[3][0]+Amat[3][3] ) - 0.5*Svec[1];

}

void Element_CPE4F::shape_functions(int nintpts, MatrixD2 & shp, MatrixD3 & dshp, VectorD & detj, MatrixD2 & x, int shp_flag )
{
	assert(nintpts == 1 || nintpts == 4);
	
	MatrixD2 s(boost::extents[nintpts][2]);
	MatrixD2 dnds(boost::extents[nen][ndm]); // temp
	MatrixD2 xjac(boost::extents[ndm][ndm]); // temp
	MatrixD2 dxds(boost::extents[ndm][ndm]); // temp
	double xsj;           		             // temp
	
	if ( nintpts == 1 ) {
		
		s[0][0] = 0.0;
		s[0][1] = 0.0;	
	
	} else if ( nintpts == 4 )  {
		
		double pt = sqrt(1.0/3.0);
		s[0][0] = -pt;
		s[0][1] = -pt;
		
		s[1][0] =  pt;
		s[1][1] = -pt;
		
		s[2][0] =  pt;
		s[2][1] =  pt;
		
		s[3][0] = -pt;
		s[3][1] =  pt;
		
	} else {
	
		cout << "compute_shape CPE4F: shouldn't be here\n";
	
	}

	std::fill(  shp.origin(),  shp.origin() +  shp.num_elements(), 0.0 );
	std::fill( dshp.origin(), dshp.origin() + dshp.num_elements(), 0.0 );
	std::fill( detj.origin(), detj.origin() + detj.num_elements(), 0.0 );
		
	for (int lquad = 0; lquad < nintpts; lquad++)
	{
		//cout << "gauss coord = (" << s[lquad][0] << "," << s[lquad][1] << endl;
		// ===================================
		// fill in shape functions N1,N2,N3,N4
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
				for ( int a = 0; a < nen; ++a ) {
					dxds[i][j] += dnds[a][i]*x[a][j];
				}
			}
		}
		
		// ==============================
		// Compute determinant of Jacobian matrix
		// ==============================
		xsj = dxds[0][0]*dxds[1][1]-dxds[1][0]*dxds[0][1];
		
		if ( xsj <=1.e-20 )
		{
			cout << "=============xsj input coordinates are==============\n";
			for ( int a = 0; a < nen; ++a ) {
				cout << x[a][0] << "\t" << x[a][1] << endl;
			}	
		}
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

