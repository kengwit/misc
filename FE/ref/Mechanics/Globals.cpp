#include "Globals.h"

namespace fem 
{
unsigned int formulation = fem::updated_lagrangian; // only in C++11
int nodes;                  // current number of nodes          (numnp)   
int elements;               // current number of elements       (numel)   
int nodes_element;          // number of nodes/element          (nen)     
//int strain_dimension;       // dimension of stress/strain       (nsdm)    
int quadrature_points;      // number of quadrature pts/element (nquad) 
int internal_dimension;     // number of internal variables  
//int shape_dimension;        // dimension of shape fn array      (nshpdm)  
int spatial_dimension;      // spatial dimension                (ndm)     
int dof_node;               // number of dof per node           (ndf)     
int materials; 							  // number of materials in database

MatrixI2 connectivity;    // element connectivity       (ix)    
MatrixD3 internal;
MatrixD2 X;     // nodal coordinates      (x)       
MatrixD2   U;             // nodal displacements    (b)   
MatrixD2   Uold;             // nodal displacements    (b)   
MatrixD2  dU;             // incremental disp       (db)      

MatrixD2 fext;
MatrixD2 fint;
MatrixD2 frct;

MatrixI2 boundary;        // boundary conditions code     (id)
MatrixD2 U_prescribed;    
MatrixD2 fext_prescribed;


vector<Element *> element_type;
vector<Material *> material_database;      // material for each element

int    step    =1;     
int    nsteps  =1;
int    NRit    =1;
int    NRitmax =1;
double NRtol   =1.e-16;
double penalty_coeff = 1.e9;		
bool   NR_converged_flag = false;
	         
EigenSpMat  KtanEigen; 
EigenVec   ResidEigen;
EigenVec     ddUEigen;
EigenVec  ddReacEigen;
			

};