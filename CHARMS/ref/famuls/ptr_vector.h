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
#ifndef PTR_VECTOR_H
#   define PTR_VECTOR_H


#include "famuls.h"
#include <cstdio>
#include <algorithm>
#include "famexception.h"

template <typename PTR, int DIM>
class PTR_VECTOR {
  
 public: // class declarations //////////////////////////////////////
  
  /**
     Construct a pointer vector; all elements set to NULL.
   */
  PTR_VECTOR () { this->zero (); }

  /**
     Produces vector with all components copied.
  */
  PTR_VECTOR (const PTR_VECTOR<PTR,DIM>& V);
  
  ~PTR_VECTOR ();
  
  /**
     Get the size of the vector.
  */
  sizet size () { return DIM; }
  
  /**
     Access element.
  */
  PTR operator() (const unsigned int i) const;
  
  /**
   Access element as writable.
   */
  PTR& operator() (const unsigned int i);
  
 private: // object methods ///////////////////////////////////////
  
  /**
     Set to NULL.
  */
  void zero ();
  
 private: // object data //////////////////////////////////////////
  
  PTR _data[DIM];

};


template <typename PTR, int DIM>
inline
PTR PTR_VECTOR<PTR,DIM>::operator() (const unsigned int i) const
{
  return _data[i];
}


template <typename PTR, int DIM>
inline
PTR & PTR_VECTOR<PTR,DIM>::operator() (const unsigned int i) 
{
  return _data[i];
}

template <typename PTR, int DIM>
PTR_VECTOR<PTR,DIM>::~PTR_VECTOR () {}


template <typename PTR, int DIM>
void PTR_VECTOR<PTR,DIM>::zero () {
  for (int j = 0; j < DIM; j++) _data[j] = 0;
}

#endif
