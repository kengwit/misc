#include "NeoHookean.h"


void NeoHookean::initialize_state_variables(ViewD1 & svars)
{
	fill(svars.origin(),svars.origin()+svars.num_elements(),0.0);
}	

void NeoHookean::constitutive_update(unsigned int el, unsigned int gp, unsigned int sec_type, 
										VectorD & props, ViewD2 & F_t, ViewD2 & F_tau, 
										ViewD2 & Dmat, ViewD1 & Svec, 
										ViewD1 & svars, bool tangent_flag)
{
	if ( fem::formulation != fem::updated_lagrangian ) {
		
		cout << "NeoHookean material requires large-strain formulation." << endl << endl;
		cout << "Set fem::formulation = fem::updated_lagrangian in main file." << endl << endl;
		assert(fem::formulation == fem::updated_lagrangian);
	}
	Tensor2 Iden(boost::extents[3][3]);
	Tensor2 F_tauT(boost::extents[3][3]);
	Tensor2 Fe_old(boost::extents[3][3]);
	Tensor2 B(boost::extents[3][3]);
	Tensor2 T_tau(boost::extents[3][3]);
	Tensor4 SpTanMod(boost::extents[3][3][3][3]);
		
	int sv_count;
	double lambda,mu,detF;
	
	// get state variables
	if ( (fem::step == 1) && (fem::NRit == 1) ) // initialize at time 0
	{
		onem(Fe_old);
		
	} else {
		
		sv_count=0;
		for ( int i = 0; i < 3; ++i ) {
			for ( int j = 0; j < 3; ++j ) {
				Fe_old[i][j] = svars[sv_count];
				sv_count++;
			}
		}
		
	}
	
	
	// identity tensor
    onem(Iden);
	
    // Obtain material properties
    lambda = props[0];
    mu     = props[1];
	
	detF = matlib_determinant(F_tau);
	
	/*cout << "F_tau:\n";
	for ( int i = 0; i < 3; i++ ) {
		for ( int j = 0; j < 3; j++ ) {
			cout << F_tau[i][j] << " ";	
		}
		cout << endl;
	}*/
	//cin.get();
	
	assert(detF > 1.e-20);
	
	// transpose of F
	matlib_transpose(F_tau,F_tauT);
	
	/*cout << "F_tauT:\n";
	for ( int i = 0; i < 3; i++ ) {
		for ( int j = 0; j < 3; j++ ) {
			cout << F_tauT[i][j] << " ";	
		}
		cout << endl;
	}*/
	//cin.get();
	
	// compute left Cauchy-Green tensor B=F*F^T
	matlib_matmult(F_tau, F_tauT, B);
	
	// compute Cauchy stress (see lecture notes of Ortiz)
	// eq. 4.7.22c
	zero(T_tau);
	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < 3; ++j )
			T_tau[i][j] = ( mu*B[i][j] + (lambda*log(detF)-mu)*Iden[i][j] )/detF;
	
	
	
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
		// compute spatial tangent moduli (see lecture notes of Ortiz)
		// eq. 4.7.26b
		zero(SpTanMod);				
		for ( int i = 0; i < 3; ++i )
			for ( int j = 0; j < 3; ++j )
				for ( int k = 0; k < 3; ++k )
					for ( int l = 0; l < 3; ++l )
						SpTanMod[i][j][k][l] = (  lambda*Iden[i][j]*Iden[k][l] + (mu-lambda*log(detF))*(Iden[i][k]*Iden[j][l]+Iden[i][l]*Iden[j][k])  )/detF + Iden[i][k]*T_tau[j][l];
		
		// Voigt form of tangent moduli
		if ( sec_type == PLANE_STRAIN ) {
			
			// note the size of Dmat is for general unsymmetric moduli, 
			// based on D.22 pg 763 in Peric complas 
			// kl should jibe with gradient operator Gmat ordering
			Dmat[0][0] = SpTanMod[0][0][0][0]; // ij=11, kl=11
			Dmat[0][1] = SpTanMod[0][0][1][0]; // ij=11, kl=21
			Dmat[0][2] = SpTanMod[0][0][0][1]; // ij=11, kl=12
			Dmat[0][3] = SpTanMod[0][0][1][1]; // ij=11, kl=22
			
			Dmat[1][0] = SpTanMod[1][0][0][0]; // ij=21, ...
			Dmat[1][1] = SpTanMod[1][0][1][0];
			Dmat[1][2] = SpTanMod[1][0][0][1];
			Dmat[1][3] = SpTanMod[1][0][1][1];
			
			Dmat[2][0] = SpTanMod[0][1][0][0]; // ij=12, ...
			Dmat[2][1] = SpTanMod[0][1][1][0];
			Dmat[2][2] = SpTanMod[0][1][0][1];
			Dmat[2][3] = SpTanMod[0][1][1][1];
			
			Dmat[3][0] = SpTanMod[1][1][0][0]; // ij=22, ...
			Dmat[3][1] = SpTanMod[1][1][1][0];
			Dmat[3][2] = SpTanMod[1][1][0][1];
			Dmat[3][3] = SpTanMod[1][1][1][1];
			
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
			svars[sv_count] = F_tau[i][j]; // Fe_old = F_tau
			sv_count++;
		}
	}
	
		
    
}
