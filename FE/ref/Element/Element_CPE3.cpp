#include "Element_CPE3.h"

// =================================================================================
//                          2D TRI 3-NODE PLANE STRAIN ELEMENT 
// =================================================================================
Element_CPE3::Element_CPE3(unsigned int id) : Element(id)
{ 
	// cout << "Destroy CPE3 element " << el_index << endl; 
	element_formulation = STANDARD;
	section_type = PLANE_STRAIN;
	nen          = 3;                   // nodes per element
	ndm          = 2;                   // number of dimensions
	//ndm1         = ndm+1;
	ndf          = 2;                   // number of dof per node
	nquad        = 1;                   // number of quadrature points
	//nshp         = (ndm+1);             // shape dimension = shape function itself plus its spatial derivatives  = 4 per node 
	nsdm         = ndm*ndf;             // strain dimension = 2*2 = 4
	//nd           = nsdm*nsdm,           // size of constitutive tangent moduli 
	nst          = ndf*nen;             // number of dof per element
	//ns           = nst*nst,             // size of K (element stiffness matrix) 	
	nvoigt       = ndm*(ndm+1)/2;       // number of stress terms in voigt notation
	
	assert(nvoigt == 3);
	
	allocate_data();
	
	
} 

Element_CPE3::~Element_CPE3() 
{ 
	// cout << "Destroy CPE3 element " << el_index << endl; 
	clear_data();
}
 


void Element_CPE3::compute_strains(int nintpts, MatrixD3 & eps, MatrixD3 & dshp, MatrixD2 & x, MatrixD2 & u, MatrixD2 & du, MatrixD2 & q )
{
	assert(nintpts == 1);
	
	// zero out first
	std::fill( eps.origin(), eps.origin() + eps.num_elements(), 0.0 );
	
	// compute deformation gradients FiJ at time n+1
	for ( int lquad = 0; lquad < nintpts; lquad++ )
	{
		
		if ( fem::formulation == fem::small_strain ) // small strain -- compute linearized strain tensor
		{
			for(int i = 0; i < ndf; i++) {
				for(int J = 0; J < ndm; J++) {
					for(int a = 0; a < nen; a++) {
						eps[lquad][i][J]+=0.5*(u[a][i]*dshp[lquad][a][J]+u[a][J]*dshp[lquad][a][i]); 
					}
				}    
			}
			
		} else if ( fem::formulation == fem::updated_lagrangian ) { // updated lagrangian
			
			for(int i = 0; i < ndf; i++) {
				for(int J = 0; J < ndm; J++) {
					for(int a = 0; a < nen; a++) {
						eps[lquad][i][J]+=(x[a][i]+u[a][i])*dshp[lquad][a][J]; 
					}
				}    
			}
			
			eps[lquad][2][2] = 1.0; // plane strain	condition
			
		} else { 
			
			cout << "unrecognized formulation in compute strain of CPE3\n";
			exit(1);
		}
				
	}
	
}	
	
void Element_CPE3::weights(int nintpts, VectorD & w)
{
	assert(nintpts == 1);
	w[0] = 1.0;

}

void Element_CPE3::compute_Bmat(int nintpts, MatrixD3 & B, MatrixD3 & dshp)
{
	std::fill(  B.origin(),  B.origin() +  B.num_elements(), 0.0 );
	
	int skip;
		
	// see Peric complas book Eq. 4.30 pg 89 - we follow convention of this book
	for (int lquad = 0; lquad < nintpts; lquad++)
	{
		for (int a = 0; a < nen; a++)
		{
			skip = a*ndf;
			B[lquad][0][skip  ] = dshp[lquad][a][0]; // d/dx -> 11
			B[lquad][1][skip+1] = dshp[lquad][a][1]; // d/dy -> 22
			B[lquad][2][skip  ] = dshp[lquad][a][1]; // d/dy -> 12
			B[lquad][2][skip+1] = dshp[lquad][a][0]; // d/dx -> 21
		}		
	}
		
}

void Element_CPE3::compute_Gmat(int nintpts, MatrixD3 & G, MatrixD3 & dshp)
{
	std::fill(  G.origin(),  G.origin() +  G.num_elements(), 0.0 );
	
	int skip;
		
	// see Peric complas book Eq. 4.94 - we follow convention of this book
	for (int lquad = 0; lquad < nintpts; lquad++)
	{
		for (int a = 0; a < nen; a++)
		{
			skip = a*ndf;
			G[lquad][0][skip  ] = dshp[lquad][a][0]; // d/dx -> 11
			G[lquad][1][skip+1] = dshp[lquad][a][0]; // d/dx -> 21
			G[lquad][2][skip  ] = dshp[lquad][a][1]; // d/dy -> 12 
			G[lquad][3][skip+1] = dshp[lquad][a][1]; // d/dy -> 22
		}		
	}
	
}

void Element_CPE3::shape_functions(int nintpts, MatrixD2 & shp, MatrixD3 & dshp, VectorD & detj, MatrixD2 & x, int shp_flag)
{

	assert(nintpts == 1);
	
	// compute shape functions and derivatives
	// r = s = 1/3
	// t = 1-r-s = 1/3
	MatrixD2 s(boost::extents[nintpts][3]);
	MatrixD2 dnds(boost::extents[nen][ndm]); // temp
	MatrixD2 xjac(boost::extents[ndm][ndm]); // temp
	MatrixD2 dxds(boost::extents[ndm][ndm]); // temp
	double xsj;           		             // temp
	
	s[0][0] = 1./3.;
	s[0][1] = 1./3.;
	s[0][2] = 1./3.;
	
	std::fill(  shp.origin(),  shp.origin() +  shp.num_elements(), 0.0 );
	std::fill( dshp.origin(), dshp.origin() + dshp.num_elements(), 0.0 );
	std::fill( detj.origin(), detj.origin() + detj.num_elements(), 0.0 );
		
	for (int lquad = 0; lquad < nintpts; lquad++)
	{
		// ===================================
		// fill in shape functions N1,N2,N3
		// ===================================
		if ( shp_flag == true )
		{
			shp[lquad][0] = s[lquad][0]; // N1 = r
			shp[lquad][1] = s[lquad][1]; // N2 = s
			shp[lquad][2] = s[lquad][2]; // N3 = t = 1-r-s
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
		// remember that u = 1-r-s-t
		dnds[0][0] =  1.0; // dN1/dr
		dnds[0][1] =  0.0; // dN1/ds

		dnds[1][0] =  0.0; // dN2/dr
		dnds[1][1] =  1.0; // dN2/ds

		dnds[2][0] = -1.0; // dN3/dr
		dnds[2][1] = -1.0; // dN3/ds

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
		for ( int i = 0; i < ndm; ++i )
			for ( int j = 0; j < ndm; ++j )
				for ( int k = 0; k < nen; ++k )
					dxds[i][j] += dnds[k][i]*x[k][j];
		
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
		detj[lquad] = xsj/2.0; // 2.0 because triangle area = det_jacobian / 2.0
		
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

