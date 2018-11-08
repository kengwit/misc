#ifndef LINEAR_ELASTIC_H
#define LINEAR_ELASTIC_H

#include "Material.h"


 

class LinearElastic : public Material
{
public:
	LinearElastic() { 
		//cout << "LinearElastic construct" << endl;
		mat_type  = "LinearElastic";
		//
		// number of properties = 2
		// lambda = NU * YOUNGS / (1+NU) / (1-2*NU);
		// mu     = YOUNGS / 2.0 / (1+NU);
		//
		nprops  = 2; 
		nsvars  = 9; // we will store at least F_old and Fe_old as state variables at each gauss point of element
		
		allocate_data();		
		
	} 
    ~LinearElastic() { 
		//cout << "LinearElastic destroy" << endl; 		
		clear_data();
	} 
	
	void initialize_state_variables(ViewD1 & svars);
    void constitutive_update(unsigned int el, unsigned int gp, unsigned int sec_type, VectorD & props, ViewD2 & eps_t, ViewD2 & eps_tau, ViewD2 & Dmat, ViewD1 & Svec, ViewD1 & svars, bool tangent_flag); 
};


// Material utility functions
void onemD2(MatrixD2 & A);
void zeroD2(MatrixD2 & A);
void zeroD4(MatrixD4 & A);


#endif