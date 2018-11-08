#include "Mechanics.h"

void mech_loop_through_elements(int isw)
{
	if ( isw == DO_KMAT_RES ) 
	{
		// zero out global internal force vector
		std::fill( fem::fint.origin(), fem::fint.origin() + fem::fint.num_elements(), 0.0 );
		
		// zero out global tangent stiffness matrix
		fem::KtanEigen.setZero(); // only removes the nonzeros, but keep allocated memory, i.e. this does not clear memory
		fem::KtanEigen.data().squeeze(); //release as much memory as possible 
		
	}

	//cout << "mechloop1\n";
	
	Element * elData;
	MatrixI2::index_gen indI2;
	MatrixD3::index_gen indD3;
	for ( int el = 0; el < fem::elements; ++el )
	{
		elData = fem::element_type[el];
		elData->isw  = isw;
		
		if ( isw == DO_SHAPE || isw == DO_STRAIN || isw == DO_KMAT_RES )
		{
			ViewI1 conl  = fem::connectivity[ indI2[el][range()] ]; // DO_SHAPE,DO_STRAIN,DO_KMAT_RES
			elData->conl = conl.origin(); // DO_SHAPE,DO_STRAIN,DO_KMAT_RES
			
			if ( isw == DO_SHAPE || isw == DO_STRAIN )
			{
				// copy ref coords and total displacements // DO_SHAPE,DO_STRAIN
				
				if ( fem::formulation == fem::small_strain )
				{
					for ( int i = 0; i < elData->nen; ++i ) {	
						for ( int j = 0; j < elData->ndm; ++j ) {
							elData->ul [i][j] = fem::U [ conl[i] ][ j ];    // total displ
							elData->dul[i][j] = fem::dU[ conl[i] ][ j ];    // incremental displ within time step
							elData->xl [i][j] = fem::X [ conl[i] ][ j ];    // reference coordinates
							
						}
					}
					
				} else if ( fem::formulation == fem::updated_lagrangian ) {
					//cout << "elData->uoldl.num_elements() = " << elData->uoldl.num_elements() << endl;
					for ( int i = 0; i < elData->nen; ++i ) {	
						for ( int j = 0; j < elData->ndm; ++j ) {
							elData->ul [i][j]   = fem::U [ conl[i] ][ j ];    // total displ
							elData->dul[i][j]   = fem::dU[ conl[i] ][ j ];    // incremental displ within time step
							elData->xl [i][j]   = fem::X [ conl[i] ][ j ];    // reference coordinates
							elData->xCl[i][j]   = elData->xl[i][j] + elData->ul[i][j]; // deformed coordinates
							//elData->uoldl[i][j] = elData->ul[i][j] - elData->dul[i][j]; // previous converged displacements
							elData->uoldl[i][j] = fem::Uold[ conl[i] ][ j ];
							//cout << "uoldl[ " << i << "][" << j << "] = " << elData->uoldl[i][j] << endl;
							//cout << "ul[ " << i << "][" << j << "] - dul[ " << i << "][" << j << "] = " <<  elData->ul[i][j] - elData->dul[i][j] << endl;
							//cin.get();
							
						}
					}
					
				} else {
					
					cout << "Mechanics.cpp - unrecognized formulation for attach data\n";
					
				}
				
			}
			
			if ( isw == DO_KMAT_RES )
			{
				//cout << "mechloop2\n";
				
				// copy history // DO_KMAT_RES
				ViewD2 ql    = fem::internal[ indD3[el][range()][range()] ]; // DO_KMAT_RES
				//cout << "ql.num_elements() = " << ql.num_elements() << endl;
				//cout << "elData->ql.num_elements() = " << elData->ql.num_elements() << endl;
				for ( int i = 0; i < elData->nquad; ++i ) {	
					for ( int j = 0; j < elData->nsvars; ++j ) {
						elData->ql[i][j] = ql[i][j];
					}
				}
				
				//cout << "mechloop3\n";
	
			}
			
		}
		
		// call element driver 		
		ElementDriver(elData);
		
		// assemble //mech_add_local(elData); // move to Element

		
	}
	
	
	
}

void mech_implicit_solve()
{
	std::fill( fem::U.origin()   , fem::U.origin() + fem::U.num_elements(), 0.0 );
	std::fill( fem::Uold.origin(), fem::Uold.origin() + fem::Uold.num_elements(), 0.0 );
	std::fill( fem::frct.origin(), fem::frct.origin() + fem::frct.num_elements(), 0.0 );
		
	double resid_norm,resid_norm0,rel_norm;
	
	// set below in main file
	//fem::nsteps  = 100;
	//fem::NRtol   = 1.e-6;
	//fem::NRitmax = 20;
	
	// ===================================================
	//  TIME STEPPING
	// ===================================================		
	for ( fem::step = 1; fem::step <= fem::nsteps; fem::step++ )
	{
		cout << "STEP: " << fem::step << endl;
		
		fem::NRit              = 1;
		fem::NR_converged_flag = false;
		resid_norm             = 2.0*fem::NRtol;
		
		// ===================================================
		//  ITERATE
		// ===================================================
		// zero out dU
		if (fem::NRit >= 1) cout << "\t     iter\t   resid\t   rel. norm " << endl;
	
		for ( int i = 0; i < fem::nodes*fem::dof_node; ++i ) {
			fem::U.origin()[i] = fem::Uold.origin()[i];
		}
	
		std::fill( fem::dU.origin(), fem::dU.origin() + fem::dU.num_elements(), 0.0 );
		while ( ( fem::NRit <= fem::NRitmax ) && ( resid_norm > fem::NRtol ) )
		{	
			// ===================================================
			//  compute shape functions and their derivatives,  
			//  strains, tangent stiffness and internal forces
			//cout << "here1\n";
			mech_loop_through_elements(DO_SHAPE);	
			//cout << "here2\n";
			mech_loop_through_elements(DO_STRAIN);	
			//cout << "here3\n";
			mech_loop_through_elements(DO_KMAT_RES);
			//cout << "here4\n";
			// ===================================================
			//  solve for incremental displacements ddU
			double prop_fac = ( double(fem::step) / double(fem::nsteps) );
			
			for ( int i = 0; i < fem::nodes*fem::dof_node; ++i ) {
				fem::fext.origin()[i] = prop_fac*fem::fext_prescribed.origin()[i];
			}
			
			for ( int i = 0; i < fem::nodes*fem::dof_node; ++i ) {
				fem::ResidEigen.coeffRef(i) =  fem::fext.origin()[i] + fem::frct.origin()[i] - fem::fint.origin()[i];
			}
		
			// solve for ddU and ...
			if ( fem::NRit == 1 ) {
				// impose non-homogeneous dirichlet boundary conditions
				// incremental prescribed displs = (1/lstps)*bc
				mech_solve_linearized_system( 1.0 / double(fem::nsteps) );				
			} else {
				// solve system for homogeneous dirichlet boundary conditions
				mech_solve_linearized_system(0.0);
			}

			
			// ===================================================
			//  update displacements, reactions and residual
			for ( int i = 0; i < fem::nodes*fem::dof_node; ++i ) {
				
				// update displacements
				fem::dU.origin()[i] += fem::ddUEigen.coeff(i); // this is incremental displacement within this time step
				fem::U.origin()[i]  += fem::ddUEigen.coeff(i); // this is total displacements up till now
				
				// update reactions
				fem::frct.origin()[i] += fem::ddReacEigen.coeff(i);

				// calculate residual with new reactions ( now fint contains nodal forces from stresses )
				fem::ResidEigen.coeffRef(i) = fem::frct.origin()[i] + fem::fext.origin()[i] - fem::fint.origin()[i];
				
			}
			
			/*cout << "--------------------------\n";
			cout << "fem::ddU = " << endl;
			for ( int i = 0; i < fem::nodes*fem::dof_node; ++i ) {
			cout << fem::ddUEigen.coeff(i) << endl;
			}
			cout << "--------------------------\n";
			cout << "fem::fint = " << endl;
			for ( int i = 0; i < fem::nodes*fem::dof_node; ++i ) {
			cout << fem::fint.origin()[i] << endl;
			}
			cout << "--------------------------\n";
			cout << "fem::frct = " << endl;
			for ( int i = 0; i < fem::nodes*fem::dof_node; ++i ) {
			cout << fem::frct.origin()[i] << endl;
			}
			cout << "--------------------------\n";
			cout << "fem::fext = " << endl;
			for ( int i = 0; i < fem::nodes*fem::dof_node; ++i ) {
			cout << fem::fext.origin()[i] << endl;
			}*/
			
			// ===================================================
			//  check convergence
			
			// compute residual norm
			resid_norm = fem::ResidEigen.norm();
			
			if ( fem::NRit == 1 ) {
				resid_norm0 = resid_norm;
				rel_norm = 1.0;
			} else {
				rel_norm = resid_norm/resid_norm0;
			}
			cout << std::scientific;
			cout.precision(5);
			cout.width(12);
			if (fem::NRit >= 1) cout << "\t" << fem::NRit << "   \t" << resid_norm << "   \t" << rel_norm << endl;
			
			if ( (fem::NRit <= fem::NRitmax) && ( resid_norm < fem::NRtol ) ) {
				fem::NR_converged_flag = true;
				break;
			}
			
			fem::NRit++;
        
		}
		
		if ( fem::NR_converged_flag == true ) {
			
			cout << "\t convergence in " << fem::NRit << " iterations\n";
			cout << "\t fem::U = " << endl;
			cout << std::scientific;
			cout.precision(5);
			for ( int i = 0; i < fem::U.shape()[0]; ++i ) {
				cout << "\t";
				for ( int j = 0; j < fem::U.shape()[1]; ++j ) {
					cout.width(12); cout << fem::U[i][j] << " ";
				}
				cout << endl;
			}
			
			cout << "\t fem::X = " << endl;
			cout << std::scientific;
			cout.precision(5);
			for ( int i = 0; i < fem::X.shape()[0]; ++i ) {
				cout << "\t";
				for ( int j = 0; j < fem::X.shape()[1]; ++j ) {
					cout.width(12); cout << fem::X[i][j]+fem::U[i][j] << " ";
				}
				cout << endl;
			}
			
			mech_loop_through_elements(DO_COMMIT_STATE);
		
			for ( int i = 0; i < fem::nodes*fem::dof_node; ++i ) {
				fem::Uold.origin()[i] = fem::U.origin()[i];
			}
			
		} else {
			cout << "no convergence\n";  
		}
    
	} 

	
}

void mech_solve_linearized_system(double scale_disp_bc)
{
	// penalty approach for now
	
	
	// =====================================================
	// copy appliedBC
	// =====================================================
	Eigen::VectorXd ddUBc(fem::nodes*fem::dof_node);
	//cout << "scale_disp_bc: " << scale_disp_bc << endl;
	for ( int i = 0; i < fem::nodes*fem::dof_node; i++ ) {
		ddUBc.coeffRef(i) = scale_disp_bc*fem::U_prescribed.origin()[i];
		//cout << "ddUBc.coeffRef(" << i << ") = " << ddUBc.coeffRef(i) << endl;
	}
	
	
	double CLarge=fem::KtanEigen.diagonal().maxCoeff()*fem::penalty_coeff;
	
	//EigenSpMat  KtanEigenModified(fem::KtanEigen);
	//KtanEigenModified.makeCompressed();
	
	// =====================================================
	// modify 
	// =====================================================
	for ( int i = 0; i < fem::nodes*fem::dof_node; i++ ) {
		int code = fem::boundary.origin()[i];		
		if ( code == DISPL_BC )
		{
			//cout << "fixed i = " << i << endl;
			//KtanEigenModified.coeffRef(i,i) += CLarge;	
			fem::KtanEigen.coeffRef(i,i) += CLarge;	
			fem::ResidEigen.coeffRef(i)  += CLarge*ddUBc.coeff(i);
			
		}
		
	}
	
	/*cout << "ke = " << endl;
	for ( int i = 0; i < 8; i++ ) {
		for ( int j = 0; j < 8; j++ ) {
			cout.width(20);
			cout.precision(12);
			
			cout << std::right << fem::KtanEigen.coeff(i,j) << " ";
		}
		cout << endl;
	}*/
	
	/*cout << "rhs = " << endl;
	for ( int i = 0; i < 8; i++ ) {
			cout.width(20);
			cout.precision(12);
			
			cout << std::right << fem::ResidEigen.coeff(i) << endl;
	}*/
	
	// solve for ddU
	fem::KtanEigen.makeCompressed();
	
	Eigen::SparseLU<EigenSpMat> solver; // SparseLU works for unsymmetric (manual says input matrix must be compressed and column major)
	
	//solver.compute(KtanEigenModified);
	solver.compute(fem::KtanEigen);
	if(solver.info()!=Eigen::Success) {
		cout << " decomposition failed " << endl;
		return;
	}
	fem::ddUEigen = solver.solve(fem::ResidEigen);
	if(solver.info()!=Eigen::Success) {
		cout << " solving failed " << endl;
		return;
	}
	
	//cout << "fem::ddUEigen: " << endl << fem::ddUEigen << endl;
	// solve for incremental reactions
	fem::ddReacEigen.setZero();
	for ( int i = 0; i < fem::nodes*fem::dof_node; i++ ) {
		int code = fem::boundary.origin()[i];		
		if ( code == DISPL_BC ) {
			fem::ddReacEigen.coeffRef(i) = -CLarge*(fem::ddUEigen.coeff(i) - ddUBc.coeff(i));
		}		
					
		//cout << "ddREigen.coeffRef(" << i << "): " << fem::ddReacEigen.coeff(i) << endl;

	}
	
	
}


/*void mech_implicit_direct_solve_incomplete()
{
	std::fill( fem::U.origin()   , fem::U.origin() + fem::U.num_elements(), 0.0 );
	std::fill( fem::Uold.origin(), fem::Uold.origin() + fem::Uold.num_elements(), 0.0 );
	std::fill( fem::frct.origin(), fem::frct.origin() + fem::frct.num_elements(), 0.0 );
		
	double resid_norm;
	
	fem::nsteps  = 1;
	fem::NRtol   = 1.e-12;
	fem::NRitmax = 15;
	
	// ===================================================
	//  TIME STEPPING
	// ===================================================		
	for ( fem::step = 1; fem::step <= fem::nsteps; fem::step++ )
	{
		cout << "STEP: " << fem::step << endl;
		
		fem::NRit              = 1;
		fem::NR_converged_flag = false;
		resid_norm             = 2.0*fem::NRtol;
		
		// ===================================================
		//  ITERATE
		// ===================================================
		// zero out dU
		if (fem::NRit >= 1) cout << "\t     iter\t   resid " << endl;
	
		for ( int i = 0; i < fem::nodes*fem::dof_node; ++i ) {
			fem::U.origin()[i] = fem::Uold.origin()[i];
		}
	
		std::fill( fem::dU.origin(), fem::dU.origin() + fem::dU.num_elements(), 0.0 );
		while ( ( fem::NRit <= fem::NRitmax ) && ( resid_norm > fem::NRtol ) )
		{	
			cout << fem::NRit << endl;
			// ===================================================
			//  compute shape functions and their derivatives,  
			//  strains, tangent stiffness and internal forces
			mech_loop_through_elements(DO_SHAPE);	
			mech_loop_through_elements(DO_STRAIN);	
			mech_loop_through_elements(DO_KMAT_RES);
			// ===================================================
			//  solve for incremental displacements ddU
			double prop_fac = ( double(fem::step) / double(fem::nsteps) );
			
			for ( int i = 0; i < fem::nodes*fem::dof_node; ++i ) {
				fem::fext.origin()[i] = prop_fac*fem::fext_prescribed.origin()[i];
			}
			
			for ( int i = 0; i < fem::nodes*fem::dof_node; ++i ) {
				fem::ResidEigen.coeffRef(i) =  fem::fext.origin()[i] + fem::frct.origin()[i] - fem::fint.origin()[i];
				cout << "fem::ResidEigen.coeffRef(i) = " << fem::ResidEigen.coeffRef(i) << endl;
				cout << "fem::fext = " << fem::fext.origin()[i] << endl;
				cout << "fem::frct = " << fem::frct.origin()[i] << endl;
			}
		
			cout << "fem::fint = " << endl;
			for ( int i = 0; i < fem::nodes*fem::dof_node; ++i ) {
				cout << fem::fint.origin()[i] << endl;
			}
			cin.get();
			// build the 4 matrices
			EigenSpMat K11; 
			EigenSpMat K12; 
			EigenSpMat K21; 
			EigenSpMat K22; 
			
			EigenVec f1;
			EigenVec dp;
			
			vector<int> free_dofs_vec;
			vector<int> prescribed_dofs_vec;
			
			//cout << "fem::U_prescribed = " << endl;
			for ( int i = 0; i < fem::nodes*fem::dof_node; i++ )
			{
				//cout << fem::boundary.origin()[i] << endl;
				if ( fem::boundary.origin()[i] == DISPL_BC )
					prescribed_dofs_vec.push_back(i);
				else
					free_dofs_vec.push_back(i);				
			}
			
			
			
			int n1 = free_dofs_vec.size();
			int n2 = prescribed_dofs_vec.size();
			
			
			cout << "free_dofs_vec = "  << endl;
			for ( int i = 0; i < n1; i++ )
				cout << free_dofs_vec[i] << endl;
			
			cout << "prescribed_dofs_vec = " << endl;
			for ( int i = 0; i < n2; i++ )
				cout << prescribed_dofs_vec[i] << endl;
			
			K11.resize(n1,n1); K11.reserve(n1*n1);
			K12.resize(n1,n2); K12.reserve(n1*n2);
			K21.resize(n2,n1); K21.reserve(n2*n1);
			K22.resize(n2,n2); K22.reserve(n2*n2);
			f1.resize(n1);
			dp.resize(n2);
			
			
			for ( int i = 0; i < n1; i++ ) {
				int row = free_dofs_vec[i];
				for ( int j = 0; j < n1; j++ ) {
					int col = free_dofs_vec[j];
					K11.coeffRef(i,j) = fem::KtanEigen.coeff(row,col);
				}
				
			}
			
			for ( int i = 0; i < n1; i++ ) {
				for ( int j = 0; j < n2; j++ ) {
					int row = free_dofs_vec[i];
					int col = prescribed_dofs_vec[j];
					K12.coeffRef(i,j) = fem::KtanEigen.coeff(row,col);
				}
			}			
			
			for ( int i = 0; i < n2; i++ ) {
				for ( int j = 0; j < n1; j++ ) {
					int row = prescribed_dofs_vec[i];
					int col = free_dofs_vec[j];
					K21.coeffRef(i,j) = fem::KtanEigen.coeff(row,col);
				}
			}			
			
			for ( int i = 0; i < n2; i++ ) {
				int row = prescribed_dofs_vec[i];
				for ( int j = 0; j < n2; j++ ) {
					int col = prescribed_dofs_vec[j];
					K22.coeffRef(i,j) = fem::KtanEigen.coeff(row,col);
				}
				
			}	
			
			cout << "K11 = " << endl;
			for ( int i = 0; i < n1; i++ ) {
				for ( int j = 0; j < n1; j++ ) {
					cout << K11.coeff(i,j)<< " ";
				}
				cout << endl;
			}
			
			cout << "K12 = " << endl;
			for ( int i = 0; i < n1; i++ ) {
				for ( int j = 0; j < n2; j++ ) {
					cout << K12.coeff(i,j)<< " ";
				}
				cout << endl;
			}
			
			int iter1_flag = 0;
			if ( fem::NRit == 1 ) iter1_flag = 1;
			
			cout << "here1\n";
			// copy boundary conditions
			for ( int i = 0; i < n2; i++ )
			{
				int row = prescribed_dofs_vec[i];
				dp.coeffRef(i) = fem::U_prescribed.origin()[row];
				//cout << "dp = " << dp.coeff(i) << endl;
			}
			
			f1.setZero();
			for ( int i = 0; i < n1; i++ )
			{
				int row = free_dofs_vec[i];
				f1.coeffRef(i) = fem::ResidEigen.coeff(row);
				for ( int j = 0; j < n2; j++ ) {
					f1.coeffRef(i) -= (double)(iter1_flag)*K12.coeffRef(i,j)*dp.coeff(j);
				}
			}
			
			cout << "f1 = " << endl;
			for ( int i = 0; i < n1; i++ )
				cout << f1.coeff(i) << endl;
			
			K11.makeCompressed();
			K12.makeCompressed();
			K21.makeCompressed();
			K22.makeCompressed();
			
			// solve for s
			Eigen::SparseLU<EigenSpMat> solver; // SparseLU works for unsymmetric (manual says input matrix must be compressed and column major)
	
			solver.compute(K11);
			if(solver.info()!=Eigen::Success) {
				cout << " decomposition failed " << endl;
				return;
			}
			EigenVec svec = solver.solve(f1);
	
			if ( fem::NRit > 1 ) dp.setZero();
				
			// assign to displacement vector
			fill(fem::dU.origin(),fem::dU.origin()+fem::dU.num_elements(),0.0);
			
			for ( int i = 0; i < n1; i++ )
			{
				int row = free_dofs_vec[i];
				fem::dU.origin()[row] = svec.coeff(i);
			}
			
			for ( int i = 0; i < n2; i++ )
			{
				int row = prescribed_dofs_vec[i];
				fem::dU.origin()[row] = dp.coeff(i);
			}
			
			cout << "fem::dU = " << endl;
			for ( int i = 0; i < fem::dU.num_elements(); i++ )
				cout << fem::dU.origin()[i] << endl;
			
			cin.get();
			fem::NRit++;
		}
		
		if ( fem::NR_converged_flag == true ) {
			
			
		} else {
			cout << "no convergence\n";  
		}
    
	} 

	
}
*/