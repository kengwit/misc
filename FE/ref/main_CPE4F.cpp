#include "FEM.h"


void setup_model();
void clear_model();

int main()
{
	setup_model();
	
	mech_implicit_solve();
	
	clear_model();
	
	return 0;
	
}

void setup_model()
{
	fem::nsteps  = 100;
	fem::NRtol   = 1.e-8;
	fem::NRitmax = 20;
	
	
	fem::formulation = fem::updated_lagrangian;
	fem::spatial_dimension  = 2;
	fem::nodes              = 4;
	fem::elements           = 1;
	fem::nodes_element      = 4;
	fem::quadrature_points  = 4;
	fem::dof_node           = 2;
	fem::internal_dimension = 9; // Fe_old
	
	int nen,nquad,nelems,ndf,nshp,ndm,nodes,nsvars;
	
	ndm    = fem::spatial_dimension;
	nelems = fem::elements;
	ndf    = fem::dof_node;
	nen    = fem::nodes_element;
	nquad  = fem::quadrature_points;
	nodes  = fem::nodes;
	nsvars = fem::internal_dimension;
	
	
	fem::connectivity.resize(boost::extents[nelems][nen]);
	fem::X           .resize(boost::extents[nodes][ndm]);
	fem::internal    .resize(boost::extents[nelems][nquad][nsvars]);
	fem::U           .resize(boost::extents[nodes][ndf]);
	fem::Uold        .resize(boost::extents[nodes][ndf]);
	fem::dU          .resize(boost::extents[nodes][ndf]);
	
	fem::fext    .resize(boost::extents[nodes][ndf]);
	fem::fint    .resize(boost::extents[nodes][ndf]);
	fem::frct   .resize(boost::extents[nodes][ndf]);
	
	// ======================================================================================
	// boundary conditions related
	// ======================================================================================
	fem::boundary.resize(boost::extents[nodes][ndf]);
	fem::U_prescribed.resize(boost::extents[nodes][ndf]);
	fem::fext_prescribed.resize(boost::extents[nodes][ndf]);
	
	// ======================================================================================
	// solver related
	// ======================================================================================
	fem::KtanEigen.resize(nodes*ndf,nodes*ndf); fem::KtanEigen.reserve(nodes*ndf*nodes*ndf);
	fem::ResidEigen.resize(nodes*ndf);
	fem::ddUEigen.resize(nodes*ndf);
	fem::ddReacEigen.resize(nodes*ndf);
	
	// coordinates
	
	double * X = fem::X.data();
	X[0] = 0.0;
	X[1] = 0.0;
	
	X[2] = 1.0;
	X[3] = 0.0;
	
	X[4] = 1.0;
	X[5] = 1.0;
	
	X[6] = 0.0;
	X[7] = 1.0;
	
	
	// connectivity
	int * conn = fem::connectivity.data();
	conn[0] = 0;
	conn[1] = 1;
	conn[2] = 2;
	conn[3] = 3;
	
	// element type for each element
	fem::element_type.resize(nelems);
	for (int i = 0; i < nelems; ++i ) {
		fem::element_type[i] = Element::Factory(i,CPE4F); 
	}
	
	// set up two types of materials 
	fem::materials = 1;	
	fem::material_database.resize(fem::materials);
	fem::material_database[0] = Material::Factory(NEOHOOKEAN); 
	fem::material_database[0]->props[0] = 1.e4; // lambda
	fem::material_database[0]->props[1] = 1.0; // mu
	
	// assign material ids to elements
	for (int i = 0; i < nelems; ++i ) {
		fem::element_type[i]->matdb_index = 0;
		fem::element_type[i]->allocate_state_variables(); // allocate for internal variables
	}
	
	// boundary conditions 
	MatrixI2 & BC = fem::boundary;
	// node 1
	BC[0][0] = DISPL_BC;
	BC[0][1] = DISPL_BC;
	
	// node 3
	BC[3][0] = DISPL_BC;
	
	
	MatrixD2 & USet = fem::U_prescribed;
	USet[0][0] = 0.0;
	USet[0][1] = 0.0;
	USet[3][0] = 0.0;
	
	
	// initial external load
	MatrixD2 & F0 = fem::fext_prescribed;
	F0[1][0] = 1.0;
	F0[2][0] = 1.0;
	
}

void clear_model()
{
	for (int i = 0; i < fem::elements; ++i ) {
		delete fem::element_type[i]; fem::element_type[i] = NULL;
	}
	
	for (int i = 0; i < fem::materials; ++i ) {
		delete fem::material_database[i]; fem::material_database[i] = NULL;
	}
	
}

/*cout << "do shape\n";
	mech_loop_through_elements(DO_SHAPE);	
	
	Element * elData;
	for ( int el = 0; el < fem::elements; ++el )
	{
		elData = fem::element_type[el];
		
		// shape func's ref
		MatrixD2 & N = elData->shpl;
		MatrixD3 & DN = elData->dshpl;
		for ( int lquad = 0; lquad < elData->nquad; ++lquad )
			cout << N[lquad][0] << ", " << N[lquad][1] << ", " << N[lquad][2] << endl;
		
		for ( int lquad = 0; lquad < elData->nquad; ++lquad ) {
			for ( int i = 0; i < elData->nen; ++i )
			cout << DN[lquad][i][0] << ", " << DN[lquad][i][1] << endl;
			
		}
	}
	cout << "do strain\n";
	mech_loop_through_elements(DO_STRAIN);	
	
	cout << "do kmat res\n";
	mech_loop_through_elements(DO_KMAT_RES);

	Element * elData;
	std::cout << std::scientific;
	std::cout.precision(4);
	for ( int el = 0; el < fem::elements; ++el )
	{
		elData = fem::element_type[el];
		
		MatrixD2 & ke = elData->ke;	
		cout << "KE of Element " << el+1 << endl;
		for ( int i = 0; i < ke.shape()[0]; i++ ) {
			for ( int j = 0; j < ke.shape()[1]; j++ ) {
				std::cout.width(12); 
				std::cout << std::right << ke[i][j] << " ";
			}
			cout << endl;
		}
			
	}
	
	mech_solve_system();*/
	
	