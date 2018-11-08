#ifndef ELEMENT_CPE4_H
#define ELEMENT_CPE4_H

#include "Element.h"

// =================================================================================
//                          2D 4-NODE PLANE STRAIN ELEMENT 
// =================================================================================
class Element_CPE4 : public Element
{
public:
	Element_CPE4(unsigned int id);
    ~Element_CPE4();
 
	void compute_strains(int nintpts, MatrixD3 & F, MatrixD3 & dshp, MatrixD2 & x, MatrixD2 & u, MatrixD2 & du, MatrixD2 & q);
    void shape_functions(int nintpts, MatrixD2 & shp, MatrixD3 & dshp, VectorD & detj, MatrixD2 & x, int shp_flag); 
    void weights(int nintpts, VectorD & w);
	void compute_Bmat(int nintpts, MatrixD3 & B, MatrixD3 & dshp);
	void compute_Gmat(int nintpts, MatrixD3 & G, MatrixD3 & dshp);
	
	
};




#endif