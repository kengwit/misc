#include "Element.h"
#include "Element_CPE3.h"
#include "Element_CPE4.h"
#include "Element_CPE4F.h"

// =================================================================================
//                                ELEMENT DRIVER
// =================================================================================

void ElementDriver(Element * elData)
{
	
	switch(elData->isw) {
		
		case DO_SHAPE:
			
			if ( fem::formulation == fem::small_strain ) {
			
				elData->shape_functions( elData->nquad, elData->shpl, elData->dshpl, elData->detjl, elData->xl, true ); // element specific
				elData->compute_Bmat( elData->nquad, elData->Bmatl, elData->dshpl );
				elData->compute_Gmat( elData->nquad, elData->Gmatl, elData->dshpl );
				
			} else if ( fem::formulation == fem::updated_lagrangian ) {
			
				//cout << "dshp in reference -------------\n";
				elData->shape_functions( elData->nquad, elData->shpl, elData->dshpl, elData->detjl, elData->xl, true ); // reference
				//cout << "dshp in current -------------\n";
				elData->shape_functions( elData->nquad, elData->shpl, elData->dshpCl, elData->detjCl, elData->xCl, false ); // current (C)
				
				//cout << "dshp0 in reference -------------\n";
				elData->shape_functions( 1, elData->shp0l, elData->dshp0l, elData->detj0l, elData->xl, true );     // centroid (0), reference
				
				//cout << "dshp0 in current -------------\n";
				elData->shape_functions( 1, elData->shp0l, elData->dshp0Cl, elData->detj0Cl, elData->xCl, false ); // centroid (0), current (C)
				
				elData->compute_Bmat( elData->nquad, elData->Bmatl, elData->dshpCl );
				elData->compute_Gmat( elData->nquad, elData->Gmatl, elData->dshpCl );
				
				elData->compute_Gmat( 1, elData->G0matl, elData->dshp0Cl );
				
			
			} else {
			
				cout << "unrecognized formulation in ElementDriver ... exiting\n";
				exit(1);
			
			}
			
			elData->weights(elData->nquad, elData->wl); // element specific
			
			break;
		
		case DO_STRAIN:
			
			if ( fem::formulation == fem::small_strain ) {
				
				elData->compute_strains( elData->nquad, elData->Fl, elData->dshpl, elData->xl, elData->ul, elData->dul, elData->ql ); // element specific
				
			} else if ( fem::formulation == fem::updated_lagrangian ) {
			
				// current & previous deformation gradients at normal gauss points
				elData->compute_strains( elData->nquad, elData->Fl, elData->dshpl, elData->xl, elData->ul, elData->dul, elData->ql ); // element specific
				elData->compute_strains( elData->nquad, elData->Fl_t, elData->dshpl, elData->xl, elData->uoldl, elData->dul, elData->ql ); // element specific
			
				// current & previous deformation gradients at centroid
				elData->compute_strains( 1, elData->Fcl, elData->dshp0l, elData->xl, elData->ul, elData->dul, elData->ql ); // element specific
				elData->compute_strains( 1, elData->Fcl_t, elData->dshp0l, elData->xl, elData->uoldl, elData->dul, elData->ql ); // element specific			
				
			} else {
			
				cout << "unrecognized formulation in ElementDriver ... exiting\n";
				exit(1);
			
			}
			
			break;
		
		case DO_KMAT_RES:
			//cout << "comp_stress\n";
			elData->compute_stress(true,true);
			//cout << "assemble\n";
			elData->assemble(true,true);
			//cout << "done assemble\n";
			break;		
		
		case DO_COMMIT_STATE:
			elData->commit_state();
			break;
			
		default:
			cout << "unrecognized task";
			
	} // end switch
	

}

// =========================== BASE ELEMENT ==================================
Element * Element::Factory(unsigned int id, int choice)
{
	switch (choice) {
	case CPE3:
		return new Element_CPE3(id);
		break;
	
	case CPE4:
		return new Element_CPE4(id);
		break;
		
	case CPE4F:
		return new Element_CPE4F(id);
		break;
	
	default:
		cout << "Undefined element type";
		return NULL;
	}

}

void Element::allocate_state_variables()
{
	nsvars = fem::material_database[matdb_index]->nsvars;
	ql.resize(boost::extents[nquad][nsvars]);  // internal variables			

	MatrixD2::index_gen indD2;
	for ( int lquad = 0; lquad < nquad; ++lquad )
	{
		ViewD1 svars = ql[ indD2[lquad][range()] ];
		fem::material_database[matdb_index]->initialize_state_variables(svars);
	}
}

void Element::allocate_data() 
{
	// ------------------------------------------------------
	// ref and current coordinates
	// ------------------------------------------------------
	xl.resize( boost::extents[nen][ndm]); 
	xCl.resize( boost::extents[nen][ndm]);     
	
	// ------------------------------------------------------
	// for F-bar
	// shape functions & derivatives at centroid 
	// ------------------------------------------------------
	Fcl.resize(boost::extents[nquad][3][3]); // full 3x3
	Fcl_t.resize(boost::extents[nquad][3][3]); // full 3x3
	
	shp0l.resize(boost::extents[1][nen]);
		
	dshp0l.resize(boost::extents[1][nen][ndm]);
	detj0l.resize(boost::extents[1]);
		
	dshp0Cl.resize(boost::extents[1][nen][ndm]);
	detj0Cl.resize(boost::extents[1]);
	
	// for F-bar
	G0matl.resize(boost::extents[1][nsdm][nen*ndf]);
	Qmatl.resize(boost::extents[nquad][nsdm][nsdm]);
	
	// ------------------------------------------------------
	// shape functions & derivatives at normal gauss points
	// ------------------------------------------------------
	shpl.resize(boost::extents[nquad][nen]);
		
	dshpl.resize(boost::extents[nquad][nen][ndm]);
	detjl.resize(boost::extents[nquad]);
		
	dshpCl.resize(boost::extents[nquad][nen][ndm]);
	detjCl.resize(boost::extents[nquad]);
	
	wl.resize(boost::extents[nquad]);
	
	// ------------------------------------------------------
	// displacements and strains
	// ------------------------------------------------------
	dul.resize( boost::extents[nen][ndm]);  		
	ul.resize( boost::extents[nen][ndm]);      		
    uoldl.resize( boost::extents[nen][ndm]);      		
    Fl.resize(boost::extents[nquad][3][3]); // full 3x3
	Fl_t.resize(boost::extents[nquad][3][3]); // full 3x3
	 
	// ------------------------------------------------------
	// element stiffness matrix and internal force vector
	// ------------------------------------------------------
	ke.resize(boost::extents[nst][nst]);      
	fintl.resize(boost::extents[nen][ndf]);      
	
	// ------------------------------------------------------------------
	// gradient matrices, tangent moduli and stress vector (voigt-style)
	// ------------------------------------------------------------------
	Bmatl.resize(boost::extents[nquad][nvoigt][nen*ndf]); // rows = 3 for 2d (plane-strain/stress), 6 for 3d
	Gmatl.resize(boost::extents[nquad][nsdm][nen*ndf]);
	Amatl.resize(boost::extents[nquad][nsdm][nsdm]);      // note the size of Dmat is for general unsymmetric moduli, based on D.22 pg 763 in Peric complas       
	Svecl.resize(boost::extents[nquad][nvoigt]);          // symmetric cauchy stress
	
		
}

void Element::clear_data()
{
	
}
	
void Element::compute_stress(bool tangent_flag, bool res_flag)
{
	// zero out element internal force 
	std::fill( fintl.origin(), fintl.origin() + fintl.num_elements(), 0.0 );
	
	// zero out element stiffness matrix
	std::fill( ke.origin(), ke.origin() + ke.num_elements(), 0.0 );
	
	MatrixD2::index_gen indD2;
	MatrixD3::index_gen indD3;
	MatrixD5::index_gen indD5;
	
	for ( int lquad = 0; lquad < nquad; ++lquad )
	{
		// strain tensor, internal variables
		ViewD2  F_tau  =    Fl[ indD3[lquad][range()][range()] ];
		ViewD2  F_t    =    Fl_t[ indD3[lquad][range()][range()] ];
		ViewD1  svars  =    ql[ indD2[lquad][range()] ];
		
		// for voigt style tangent moduli and stress
		ViewD2   Amat  = Amatl[ indD3[lquad][range()][range()] ]; // matrix form of tangent moduli
		ViewD1   Svec  = Svecl[ indD2[lquad][range()] ];
		
		// process kinematics for F-bar
		if ( element_formulation == FBAR ) {
			ViewD2  Fc_tau  =    Fcl[ indD3[0][range()][range()] ];
			ViewD2  Fc_t    =    Fcl_t[ indD3[0][range()][range()] ];			
			process_FBarKinematics(Fc_t, Fc_tau, F_t, F_tau);
		}
		
	    // constitutive update - return stress and tangent material moduli	
		MaterialDriver(matdb_index, el_index, section_type, lquad, F_t, F_tau, Amat, Svec, svars, tangent_flag);	
		
		if ( element_formulation == FBAR ) {
			ViewD2 Qmat = Qmatl[ indD3[lquad][range()][range()] ]; // matrix form of tangent moduli
			compute_Qmat(Qmat,Amat,Svec);
		}
		
		if ( res_flag == true )
		{
			
			if ( fem::formulation == fem::small_strain ) {
				
				// Swansea style
				for ( int a = 0; a < nen; a++ ) // loop over nodes
				{
					for ( int i = 0; i < ndf; ++i ) { // loop over degrees of freedom
						for ( int j = 0; j < nvoigt; ++j ) { // loop over stress terms
							fintl[a][i] +=  Bmatl[lquad][j][ndf*a+i]*Svecl[lquad][j]*detjl[lquad]*wl[lquad];
						}
					} 
					
				} // end loop over nodes 
					
				
			} else if ( fem::formulation == fem::updated_lagrangian ) {
				
				for ( int a = 0; a < nen; a++ ) // loop over nodes
				{
					for ( int i = 0; i < ndf; ++i ) { // loop over degrees of freedom
						for ( int j = 0; j < nvoigt; ++j ) { // loop over stress terms
							fintl[a][i] +=  Bmatl[lquad][j][ndf*a+i]*Svecl[lquad][j]*detjCl[lquad]*wl[lquad];
						}
					} 
					
				} // end loop over nodes 
			
			} else {
				cout << "unrecognized formulation - internal force assembly (tensor style) ... exiting\n";
				exit(1);
			} 
			
		
			
			
		} // end if do residual
		
		
		if ( tangent_flag == true )
		{
			
			if ( fem::formulation == fem::small_strain )
			{
				// Swansea style
				for (int a=0; a<nen; a++) { // loop over nodes
					for(int i=0; i<ndf; i++){ // loop over dofs of node
						for(int b=0; b<nen; b++) { // loop over nodes
							for(int j=0; j<ndf; j++) { // loop over dofs of node					
								for(int I=0; I<nsdm; I++) { // loop over nsdm == ndm*ndf (tangent)
									for(int J=0; J<nsdm; J++) { // loop over nsdm == ndm*ndf (tangent)
										ke[a*ndf+i][b*ndf+j]+=Gmatl[lquad][I][ndf*a+i]*Amatl[lquad][I][J]*Gmatl[lquad][J][ndf*b+j]*detjl[lquad]*wl[lquad];
									}
								}
							}
						}
					}
				}
					
			} else if ( fem::formulation == fem::updated_lagrangian ) {
				
				for (int a=0; a<nen; a++) { // loop over nodes
					for(int i=0; i<ndf; i++){ // loop over dofs of node
						for(int b=0; b<nen; b++) { // loop over nodes
							for(int j=0; j<ndf; j++) { // loop over dofs of node					
								for(int I=0; I<nsdm; I++) { // loop over nsdm == ndm*ndf (tangent)
									for(int J=0; J<nsdm; J++) { // loop over nsdm == ndm*ndf (tangent)
										ke[a*ndf+i][b*ndf+j]+=Gmatl[lquad][I][ndf*a+i]*Amatl[lquad][I][J]*Gmatl[lquad][J][ndf*b+j]*detjCl[lquad]*wl[lquad];
									}
								}
							}
						}
					}
				}

				// add F-bar term if required	
				if ( element_formulation == FBAR ) {
					/*cout << "Qmat at integration point " << lquad << endl;
					cout << "Qmat\n";
					for ( int i = 0; i < 4; i++ ) {
						for ( int j = 0; j < 4; j++ ) {
							cout << Qmat[i][j] << " ";
						}
						cout << endl;
					}
					cout << "Qmatl\n";
					for ( int i = 0; i < 4; i++ ) {
						for ( int j = 0; j < 4; j++ ) {
							cout << Qmatl[lquad][i][j] << " ";
						}
						cout << endl;
					}*/
					//cin.get();
					for (int a=0; a<nen; a++) { // loop over nodes
						for(int i=0; i<ndf; i++){ // loop over dofs of node
							for(int b=0; b<nen; b++) { // loop over nodes
								for(int j=0; j<ndf; j++) { // loop over dofs of node					
									for(int I=0; I<nsdm; I++) { // loop over nsdm == ndm*ndf (tangent)
										for(int J=0; J<nsdm; J++) { // loop over nsdm == ndm*ndf (tangent)
											ke[a*ndf+i][b*ndf+j]+=Gmatl[lquad][I][ndf*a+i]*Qmatl[lquad][I][J]*( G0matl[0][J][ndf*b+j]-Gmatl[lquad][J][ndf*b+j] )*detjCl[lquad]*wl[lquad];
										}
									}
								}
							}
						}
					}
				
				}
		
			
			} else {
				cout << "unrecognized formulation - stiffness matrix assembly (tensor style) ... exiting\n";
				exit(1);
			} 
			
			
		} // end if do tangent
			
	} // end loop over quadrature points
	
	
}

void Element::assemble(bool tangent_flag, bool res_flag)
{
	/*cout << "ke = " << endl;
	for ( int i = 0; i < 8; i++ ) {
		for ( int j = 0; j < 8; j++ ) {
			cout.width(20);
			cout.precision(12);
			
			cout << std::right << ke[i][j] << " ";
		}
		cout << endl;
	}*/
	
	//cout << "resid\n";
	if ( res_flag == true )
	{			
		for ( int a = 0; a < nen; ++a ) {
			for ( int i = 0; i < ndf; ++i ) {
				fem::fint[ conl[a] ][i] += fintl[a][i];
			}
		}
	}
	
	//cout << "tang\n";
	
	if ( tangent_flag == true )
	{
		//cout << "ke.num_elements() = " << ke.num_elements() << endl;
		//cout << "conl = " <<endl;
		//for ( int i = 0; i < nen; ++i )
		//	cout << conl[i] << endl;
		
		int global_row,global_col;
	
		//cout << "assemble ke\n";
		for (int a=0; a < nen; a++) { // loop over nodes
			for(int i=0; i < ndf; i++){ // loop over dofs of node
				for(int b=0; b < nen; b++) { // loop over nodes
					for(int j=0; j < ndf; j++) { // loop over dofs of node
						
						global_row = conl[a]*ndf + i;						
						global_col = conl[b]*ndf + j;
						//cout << "conl["<<a<<"] = " << conl[a] << endl;
						//cout << "conl["<<b<<"] = " << conl[b] << endl;
						
						//cout << "a*ndf+i = " << a*ndf+i << ", b*ndf+j = " << b*ndf+j << endl;
						//cout << "ke["<<a*ndf+i<<"]["<<b*ndf+j<<"] = " <<ke[a*ndf+i][b*ndf+j] << endl;
						//cout << "fem::KtanEigen.coeff("<<global_row<<","<<global_col<<") = " << fem::KtanEigen.coeff(global_row,global_col)  <<endl; 
						fem::KtanEigen.coeffRef(global_row,global_col) += ke[a*ndf+i][b*ndf+j];
						//cout << "here done\n";
						//fem::KtanEigen.coeffRef(global_row,global_col) += 1.0;
					}
				}
			}
		}
		//cout << "done assemble KMAT_RES\n";

	}	
	//cout << "tang done\n";
	
}

void Element::commit_state()
{
	for ( int i = 0; i < nquad; ++i )
		for ( int j = 0; j < nsvars; ++j )
			fem::internal[el_index][i][j] = ql[i][j];	
}

