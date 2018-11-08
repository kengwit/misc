#ifndef ELEMENT_H
#define ELEMENT_H

#include "Globals.h"

// =======================================================================
// task enumeration
enum TASK
{ 
	DO_STRAIN, DO_SHAPE, DO_KMAT_RES, DO_COMMIT_STATE
};

enum ELEM_TYPE
{ 
	CPE3,CPE4,CPE4F
};

enum SECTION_TYPE
{
	PLANE_STRAIN,SOLID
};

enum ELEMENT_FORMULATION
{
	STANDARD,FBAR
};

// =======================================================================
class Element
{
	
public:
	
	unsigned int element_formulation;
	unsigned int isw;  // task switch
	unsigned int el_index; // element index (we differentiate between element "index" and element "label")
	unsigned int matdb_index;
	unsigned int section_type;
	int nen;
	int ndm;
	int ndf;
	int nquad;
	int nsdm;   // strain/stress dimension
	int nst;    // = nen*ndf;    
    int nvoigt; // number of stress terms in voigt notation
	int nsvars;
	
	MatrixD2    xl;     // element nodal reference coordinates
	MatrixD2   xCl;     // element nodal current coordinates
	MatrixD2   dul;     // element nodal incremental displacements
	MatrixD2    ul;     // element nodal displacements
	MatrixD2 uoldl;     // element nodal previous (converged) displacements
	MatrixD3    Fl;     // strains/def gradient at each integration point time n+1
	MatrixD3  Fl_t;     // strains/def gradient at each integration point time n
	MatrixD2    ql;     // internal variables
	
	MatrixD2 fintl;      // element nodal residual vector
	MatrixD2 ke;      // element stiffness matrix
	
	// ------------------------------------------------------
	// shape functions & derivatives at centroid (for F-bar)
	// ------------------------------------------------------
	MatrixD3    Fcl;     // strains/def gradient at each integration point time n+1
	MatrixD3  Fcl_t;     // strains/def gradient at each integration point time n
	
	MatrixD2   shp0l;    // only 1 gauss point (at centroid), not used for plane strain (used in G0mat for axisymm)
	
	MatrixD3  dshp0l;    // ref config, used to calculate deformation gradient at centroid
	VectorD   detj0l;    // ref config, not used
	
	MatrixD3 dshp0Cl;    // current config, used to calculate G0mat
	VectorD  detj0Cl;    // current config, not used
	
	MatrixD3   G0matl;   // centroid gradient operator in matrix form
	MatrixD3   Qmatl;    // F-bar Q matrix
	
	// ------------------------------------------------------
	// shape functions & derivatives at normal gauss points
	// ------------------------------------------------------
	MatrixD2    shpl;   // shape functions
	
	MatrixD3   dshpl;   // ref config, used to calculate deformation gradient
	VectorD    detjl;   // ref config, not used
	
	MatrixD3  dshpCl;   // current config, used for Bmat, Gmat
	VectorD   detjCl;   // current config, used for Bmat, Gmat  	
	
	VectorD       wl;   // weights
	
	MatrixD3   Bmatl;   // B-mat (plane-strain/stress)
	MatrixD3   Gmatl;   // Gradient operator in matrix form
	MatrixD3   Amatl;   // Tangent moduli in matrix form
	MatrixD2   Svecl;   // Symmetric Cauchy stress
	
	
	int * conl;  // element connectivity - attach to conl.origin() outside
	
    Element(unsigned int e) {
		el_index = e;	
		//cout << "Construct Base element " << el_index << endl; 
	} 
	
	virtual ~Element() { 
		//cout << "Destroy Base element " << el_index << endl; 
	} 
 
	// Factory method
	static Element * Factory(unsigned int id, int choice);
	
	void allocate_state_variables();
	
	//void compute_strains(); // now as virtual function so it is element specific
	
	void compute_stress(bool tangent_flag, bool res_flag);
	
	void assemble(bool tangent_flag, bool res_flag);
	
	void commit_state();
	
	virtual void compute_strains( int nintpts, MatrixD3 & eps, MatrixD3 & dshp, MatrixD2 & x, MatrixD2 & u, MatrixD2 & du, MatrixD2 & q )=0; // =0 to make abstract class (no instant of this class can be invoked)
    virtual void shape_functions( int nintpts, MatrixD2 & shp, MatrixD3 & dshp, VectorD & detj, MatrixD2 & x, int shp_flag )=0;//{ std::cout << "shape_functions Base" << endl; }
	virtual void weights( int nintpts, VectorD & w )=0;// { std::cout << "weights Base" << endl; }
	virtual void compute_Bmat(int nintpts, MatrixD3 & B, MatrixD3 & dshp)=0;
	virtual void compute_Gmat(int nintpts, MatrixD3 & G, MatrixD3 & dshp)=0;
	
	virtual void process_FBarKinematics(ViewD2 & Fc_t, ViewD2 & Fc_tau, ViewD2 & F_t, ViewD2 & F_tau) {}
	virtual void compute_Qmat(ViewD2 & Qmat, ViewD2 & Amat, ViewD1 & Svec) {}
	
protected:

	// ===================================================================
	//     data for connectivity and element nodal coordinates
	void allocate_data();
	
	void clear_data(); 
	// ===================================================================


};

void ElementDriver(Element * elData);


#endif

