#include "LinearElastic.h"


void LinearElastic::initialize_state_variables(ViewD1 & svars)
{
	fill(svars.origin(),svars.origin()+svars.num_elements(),0.0);
}	

void LinearElastic::constitutive_update(unsigned int el, unsigned int gp, unsigned int sec_type, 
										VectorD & props, ViewD2 & E_t, ViewD2 & E_tau, 
										ViewD2 & Dmat, ViewD1 & Svec, 
										ViewD1 & svars, bool tangent_flag)
{
	if ( fem::formulation != fem::small_strain ) {
		
		cout << "Linear Elastic material requires small-strain formulation." << endl << endl;
		cout << "Set fem::formulation = fem::small_strain in main file." << endl << endl;
		assert(fem::formulation == fem::small_strain);
	
	}
	
	// clear some variables first
	//fill( Dmat.origin, Dmat.origin + Dmat.num_elements(), 0.0 );
	//fill( Svec.origin, Svec.origin + Svec.num_elements(), 0.0 );
	
	
	int sv_count;
	MatrixD2 Ee_old(boost::extents[3][3]);
	
	// get state variables
	if ( (fem::step == 1) && (fem::NRit == 1) ) // initialize at time 0
	{
		fill(Ee_old.origin(),Ee_old.origin()+Ee_old.num_elements(),0.0);
		
	} else {
		
		sv_count=0;
		
		for ( int i = 0; i < 3; ++i ) {
			for ( int j = 0; j < 3; ++j ) {
				Ee_old[i][j] = svars[sv_count];
				sv_count++;
			}
		}
	}
	
	/*cout << "E_old\n";
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			cout << E_old[i][j] << " ";
		}
		cout << endl;
	}

	cout << "Ee_old\n";
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			cout << Ee_old[i][j] << " ";
		}
		cout << endl;
	}

	cout << "E_tau\n";
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			cout << E_tau[i][j] << " ";
		}
		cout << endl;
	}*/

	
	// Identity matrix
    MatrixD2 Iden(boost::extents[3][3]);
	onemD2(Iden);
	
    // Obtain relevant material properties
    double E  = props[0];
    double nu = props[1];
	
    double theta = E_tau[0][0]+E_tau[1][1]+E_tau[2][2];
    
	// Stress tensor
    MatrixD2 T_tau(boost::extents[3][3]);
	zeroD2(T_tau);
	
	double fac1 = E/(1.0+nu);
    double fac2 = theta*nu/(1.0-2.0*nu);
    for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
			T_tau[i][j] = fac1*( E_tau[i][j] + fac2*Iden[i][j] );

	// Voigt form of stress 
	if ( sec_type == PLANE_STRAIN ) {
		
		// see convention used in D.1 pg 759 in Peric complas book
		Svec[0] = T_tau[0][0]; // 11
		Svec[1] = T_tau[1][1]; // 22
		Svec[2] = T_tau[0][1]; // 12
		
	
	} else {
		cout << "unrecognized section type - need to be plane strain element\n";
		assert(sec_type == PLANE_STRAIN);
	}
	
	if ( tangent_flag == true )
	{		
		// Calculate the spatial tangent modulus
		//
		fac1 = E/(2.0*(1.0+nu));
		fac2 = E*nu/((1.0+nu)*(1.0-2.0*nu));
		MatrixD4 SpUUMod(boost::extents[3][3][3][3]);
		zeroD4(SpUUMod);
	
		for ( int i = 0; i < 3; ++i )
			for ( int j = 0; j < 3; ++j )
				for ( int k = 0; k < 3; ++k )
					for ( int l = 0; l < 3; ++l )
						SpUUMod[i][j][k][l] = fac1*( Iden[i][l]*Iden[j][k]+Iden[i][k]*Iden[j][l] ) + fac2*Iden[i][j]*Iden[k][l];
		
	
		// Voigt form of tangent moduli
		if ( sec_type == PLANE_STRAIN ) {
			
			// note the size of Dmat is for general unsymmetric moduli, 
			// based on D.22 pg 763 in Peric complas 
			// kl should jibe with gradient operator Gmat ordering
			Dmat[0][0] = SpUUMod[0][0][0][0]; // ij=11, kl=11
			Dmat[0][1] = SpUUMod[0][0][1][0]; // ij=11, kl=21
			Dmat[0][2] = SpUUMod[0][0][0][1]; // ij=11, kl=12
			Dmat[0][3] = SpUUMod[0][0][1][1]; // ij=11, kl=22
			
			Dmat[1][0] = SpUUMod[1][0][0][0]; // ij=21, ...
			Dmat[1][1] = SpUUMod[1][0][1][0];
			Dmat[1][2] = SpUUMod[1][0][0][1];
			Dmat[1][3] = SpUUMod[1][0][1][1];
			
			Dmat[2][0] = SpUUMod[0][1][0][0]; // ij=12, ...
			Dmat[2][1] = SpUUMod[0][1][1][0];
			Dmat[2][2] = SpUUMod[0][1][0][1];
			Dmat[2][3] = SpUUMod[0][1][1][1];
			
			Dmat[3][0] = SpUUMod[1][1][0][0]; // ij=22, ...
			Dmat[3][1] = SpUUMod[1][1][1][0];
			Dmat[3][2] = SpUUMod[1][1][0][1];
			Dmat[3][3] = SpUUMod[1][1][1][1];
			
		} else {
			cout << "unrecognized section type - need to be plane strain element\n";
			assert(sec_type == PLANE_STRAIN);
		}
	}	
	
	// update state variables
	// will write to global vector upon convergence
	sv_count=0;
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			svars[sv_count] = E_tau[i][j]; // Ee_old = E_tau (everything is elastic)
			sv_count++;
		}
	}
		
    
}

void onemD2(MatrixD2 & A)
{
	for ( int i = 0; i < 2; i++ ) assert ( A.shape()[i] == 3 );	
	
	zeroD2(A);
	
	for ( int i = 0; i < 3; ++i ) A[i][i] = 1.0;
	
}

void zeroD2(MatrixD2 & A)
{
	for ( int i = 0; i < 2; i++ ) assert ( A.shape()[i] == 3 );	
	
	std::fill( A.origin(), A.origin() + A.num_elements(), 0.0 );
}

void zeroD4(MatrixD4 & A)
{
	for ( int i = 0; i < 4; i++ ) assert ( A.shape()[i] == 3 );	
	
	std::fill( A.origin(), A.origin() + A.num_elements(), 0.0 );
}
