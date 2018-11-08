#ifndef ELEMENT_CPE3_H
#define ELEMENT_CPE3_H

#include "Element.h"

// =================================================================================
//                          2D 3-NODE PLANE STRAIN ELEMENT 
// =================================================================================
class Element_CPE3 : public Element
{
// Shape functions:
// [ dN1/dX dN1/dY N1
//   dN2/dX dN2/dY N2
//   dN3/dX dN3/dY N3 ]
// shape dimension = number of cols = 3
public:
	Element_CPE3(unsigned int id);
    ~Element_CPE3();
 
	void compute_strains(int nintpts, MatrixD3 & F, MatrixD3 & dshp, MatrixD2 & x, MatrixD2 & u, MatrixD2 & du, MatrixD2 & q);
    void shape_functions(int nintpts, MatrixD2 & shp, MatrixD3 & dshp, VectorD & detj, MatrixD2 & x, int shp_flag); 
    void weights(int nintpts, VectorD & w);
	void compute_Bmat(int nintpts, MatrixD3 & B, MatrixD3 & dshp);
	void compute_Gmat(int nintpts, MatrixD3 & G, MatrixD3 & dshp);
	
	
};




#endif