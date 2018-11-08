#ifndef GLOBALS_H
#define GLOBALS_H

#include <cassert>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
using namespace std;

#include "boost/multi_array.hpp"
typedef boost::multi_array<int, 1> VectorI;
typedef boost::multi_array<int, 2> MatrixI2;
typedef MatrixI2::array_view<1>::type ViewI1;

#include <Eigen/Sparse>
typedef Eigen::SparseMatrix<double> EigenSpMat; // declares a column-major sparse matrix type of double
typedef Eigen::VectorXd EigenVec;
typedef Eigen::Triplet<double> Triplet;

typedef boost::multi_array<double, 1> VectorD;
typedef boost::multi_array<double, 2> MatrixD2;
typedef boost::multi_array<double, 3> MatrixD3;
typedef boost::multi_array<double, 4> MatrixD4;
typedef boost::multi_array<double, 5> MatrixD5;
typedef MatrixD2::array_view<1>::type ViewD1;
typedef MatrixD3::array_view<2>::type ViewD2;
typedef MatrixD4::array_view<3>::type ViewD3;
typedef MatrixD5::array_view<4>::type ViewD4;


typedef boost::multi_array_types::index_range range;

//#include "Utilities.h"
#include "Element.h"
#include "Material.h"

// forward declaration
class Element;
class Material;

namespace fem 
{

enum { updated_lagrangian, total_lagrangian, small_strain };

#define DECLARE_EXTERN extern

DECLARE_EXTERN unsigned int formulation;
DECLARE_EXTERN int nodes;                  // current number of nodes          (numnp)   
DECLARE_EXTERN int elements;               // current number of elements       (numel)   
DECLARE_EXTERN int nodes_element;          // number of nodes/element          (nen)     
//DECLARE_EXTERN int strain_dimension;       // dimension of stress/strain       (nsdm)    
DECLARE_EXTERN int quadrature_points;      // number of quadrature pts/element (nquad)   
DECLARE_EXTERN int internal_dimension;     // number of internal variables per quadrature point
//DECLARE_EXTERN int shape_dimension;        // dimension of shape fn array      (nshpdm)  
DECLARE_EXTERN int spatial_dimension;      // spatial dimension                (ndm)     
DECLARE_EXTERN int dof_node;               // number of dof per node           (ndf)     
DECLARE_EXTERN int materials; 			   // number of materials in database

DECLARE_EXTERN MatrixD2  X;        // nodal coordinates      (x)       
DECLARE_EXTERN MatrixI2 connectivity;    // element connectivity       (ix)    
DECLARE_EXTERN MatrixD3     internal;        // internal variables     (q)
DECLARE_EXTERN MatrixD2            U;        // nodal displacements    (b)    
DECLARE_EXTERN MatrixD2         Uold;        // nodal displacements    (b)    
DECLARE_EXTERN MatrixD2           dU;        // incremental disp       (db)      

DECLARE_EXTERN MatrixI2 boundary;           // boundary code
DECLARE_EXTERN MatrixD2 U_prescribed;       // prescribed displ
DECLARE_EXTERN MatrixD2 fext_prescribed;       // total applied external load        
DECLARE_EXTERN MatrixD2 fext;               // current external load
DECLARE_EXTERN MatrixD2 fint;               // current internal load
DECLARE_EXTERN MatrixD2 frct;               // current reactions        


DECLARE_EXTERN vector<Element *>       element_type;
DECLARE_EXTERN vector<Material *> material_database;      // material for each element

// SOLVER RELATED
DECLARE_EXTERN int    step;     
DECLARE_EXTERN int    nsteps;
DECLARE_EXTERN int    NRit;
DECLARE_EXTERN int    NRitmax;
DECLARE_EXTERN double NRtol;
DECLARE_EXTERN double penalty_coeff;		
DECLARE_EXTERN bool   NR_converged_flag;
	         
DECLARE_EXTERN EigenSpMat  KtanEigen; // note: default is column-major like in fortran
DECLARE_EXTERN EigenVec   ResidEigen;
DECLARE_EXTERN EigenVec     ddUEigen;
DECLARE_EXTERN EigenVec  ddReacEigen;
};

#endif // GLOBALS_H


