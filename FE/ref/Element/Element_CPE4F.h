#ifndef ELEMENT_CPE4F_H
#define ELEMENT_CPE4F_H

#include "Element.h"

// =================================================================================
//                          2D 4-NODE F-BAR PLANE STRAIN ELEMENT 
// =================================================================================
class Element_CPE4F : public Element
{
public:
	Element_CPE4F(unsigned int id);
    ~Element_CPE4F();
 
	void compute_strains(int nintpts, MatrixD3 & F, MatrixD3 & dshp, MatrixD2 & x, MatrixD2 & u, MatrixD2 & du, MatrixD2 & q);
    void shape_functions(int nintpts, MatrixD2 & shp, MatrixD3 & dshp, VectorD & detj, MatrixD2 & x, int shp_flag); 
    void weights(int nintpts, VectorD & w);
	void compute_Bmat(int nintpts, MatrixD3 & B, MatrixD3 & dshp);
	void compute_Gmat(int nintpts, MatrixD3 & G, MatrixD3 & dshp);
	
	// F-bar specific methods
	void process_FBarKinematics(ViewD2 & Fc_t, ViewD2 & Fc_tau, ViewD2 & F_t, ViewD2 & F_tau);	
	void compute_Qmat(ViewD2 & Qmat, ViewD2 & Amat, ViewD1 & Svec);
	
};




#endif