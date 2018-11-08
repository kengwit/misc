#ifndef MATERIAL_H
#define MATERIAL_H

#include "Globals.h"

// =======================================================================
// material enumeration
enum MATERIALS
{
	LINEAR_ELASTIC, NEOHOOKEAN
};

class Material
{
	
public:
	string mat_type;
	int nprops;
	int nsvars;
	
	VectorD props;	
	
	Material() { 
		//cout << "Base Material construct" << endl; 
	} 
	
	virtual ~Material() { 
		//cout << "Base Material destroy" << endl; 
	} 
 
	// Factory Method
    static Material * Factory(int choice);
	
	virtual void initialize_state_variables(ViewD1 & svars)=0;
	virtual void constitutive_update(unsigned int el, unsigned int gp, unsigned int sec_type, VectorD & props, ViewD2 & eps_t, ViewD2 & eps_tau, ViewD2 & Dmat, ViewD1 & Svec, ViewD1 & svars, bool tangent_flag)=0;
	//{ std::cout << "Base constitutive update" << endl; }
	
	virtual void input() { }

protected:	
	// ===================================================================
	//     data for connectivity and element nodal coordinates
	void allocate_data() {
		//cout << "resize material props\n";
		props.resize(boost::extents[nprops]);
	}
	
	void clear_data() {
		//cout << "clear material props\n";
		
	}
	// ===================================================================
	
};
 
// Material driver
void MaterialDriver(unsigned int mat, unsigned int el, unsigned int gp, unsigned int sec_type, ViewD2 & eps_t, ViewD2 & eps, ViewD2 & Dmat, ViewD1 & Svec, ViewD1 & svars, bool tangent_flag);



#endif