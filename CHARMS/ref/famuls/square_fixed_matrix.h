/*--------------------------------------------------------------------------
--                                                                        --
--         (F)ramework for (A)daptive (MUL)tiphysics (S)imulations        --
--            Copyright (C) 2001, Petr Krysl (pkrysl@ucsd.edu).           --
--                                                                        --
--                 [Pronounce `famulus': scholar's helper]                --
--                                                                        --
--                                                                        --
--  This program is free software; you can redistribute it and/or modify  --
--  it under the terms of the GNU General Public License as published by  --
--  the Free Software Foundation; either version 2 of the License, or     --
--  (at your option) any later version.                                   --
--                                                                        --
--------------------------------------------------------------------------*/
#ifndef SQUARE_FIXED_MATRIX_H
#define SQUARE_FIXED_MATRIX_H

#include "famexception.h"
extern "C" {
#include "geP.h"
}

/**
   Square dense matrix.  Uses the Ckit library GEP.
*/
class SQUARE_FIXED_MATRIX {
  
 public: // object functions ////////////////////////////////////////

 /**
   Default constructor.  Matrix size is zero.
   */
  SQUARE_FIXED_MATRIX () {
    _gm.n = 0;
    _gm.pd = 0;
  }

  /**
  Make a square matrix with given size.
  */
  SQUARE_FIXED_MATRIX (unsigned int n) {
    CHECK (gep_initmat (&_gm, n), EXCEPTION_NULL_PTR ,;);
  }
  
  ~SQUARE_FIXED_MATRIX() { gep_dispose (&_gm); }

  /**
  Get the dimension of the matrix.
  */
  unsigned int n () const { return _gm.n; }

  /**
  Resize the matrix to have a new dimension n.  If the previous
  size of the matrix was the same, nothing happens; otherwise
  the old matrix is disposed of, and a new one is allocated.
  */
  void resize (unsigned int n) {
    if (_gm.n != 0) {
      if (_gm.n != (int) n) {
        gep_dispose (&_gm);
      } else {
        return; // same-size matrix
      }
    }
    CHECK (gep_initmat (&_gm, n), EXCEPTION_NULL_PTR ,;);
  }
      

  inline double operator() (const unsigned int i, const unsigned int j) const {
    return GEP_ELEM (_gm, i, j);
  }
  
  inline double& operator() (const unsigned int i, const unsigned int j) {
    return GEP_ELEM (_gm, i, j);
  }

  inline double elem (const unsigned int i, const unsigned int j) const {
    return GEP_ELEM (_gm, i, j);
  }

  /**
  Zero out matrix.
  */
  void zero () { gep_zero (&_gm); }
  
 private: // object data ////////////////////////////////////////////
  
  gep_dense_matrix_t _gm;
  
};

#endif
