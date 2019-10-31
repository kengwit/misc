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
#ifndef GENERAL_FIXED_VECTOR_H
#   define GENERAL_FIXED_VECTOR_H


#include <cstdio>
#include <algorithm>
#include <functional>
#include "famexception.h"

template <typename ELEM_TYPE>
class GENERAL_FIXED_VECTOR {
  
 public: // class declarations //////////////////////////////////////

  /**
     Construct an uninitialized vector. 
  */
  GENERAL_FIXED_VECTOR (int dim) : _dim(dim) {
    _data = new ELEM_TYPE[_dim];
  }
  /**
     Construct an uninitialized vector w/o any data. Need to use size() to allocate storage!
  */
  GENERAL_FIXED_VECTOR () { _dim = 0;_data = 0; }
  /**
     Use to allocate storage. Data is uninitialized!
  */
  void size(int dim) { _dim = dim; _data = new ELEM_TYPE[_dim]; }
  /**
     Construct an initialized vector.  The object given
     as argument is assigned to each element of the vector.
  */
  GENERAL_FIXED_VECTOR (int dim, ELEM_TYPE initval) : _dim(dim) {
    _data = new ELEM_TYPE[_dim];
    assign (initval);
  }
  ~GENERAL_FIXED_VECTOR ();
  
  /**
     Get the size of the vector.
  */
  sizet size () { return _dim; }
  
  /**
     Set to given value.
  */
  void assign (ELEM_TYPE v);
  
  /**
     Access element.
  */
  const ELEM_TYPE &operator() (const sizet i) const;
  
  /**
   Access element as writable.
   */
  ELEM_TYPE &operator() (const sizet i);
  
 private: // object data //////////////////////////////////////////
  
  sizet        _dim;
  ELEM_TYPE   *_data;

};


template <typename ELEM_TYPE>
inline
const ELEM_TYPE &GENERAL_FIXED_VECTOR<ELEM_TYPE>::operator() (const sizet i) const
{
  return _data[i];
}


template <typename ELEM_TYPE>
inline
ELEM_TYPE &GENERAL_FIXED_VECTOR<ELEM_TYPE>::operator() (const sizet i) 
{
  return _data[i];
}


template <typename ELEM_TYPE>
GENERAL_FIXED_VECTOR<ELEM_TYPE>::~GENERAL_FIXED_VECTOR () {
  if (_data) delete [] _data;
}


template <typename ELEM_TYPE>
void GENERAL_FIXED_VECTOR<ELEM_TYPE>::assign (ELEM_TYPE v) {
  for (sizet j = 0; j < _dim; j++) _data[j] = v;
}

#endif
