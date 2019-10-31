/*
                           OOOFS
                 Copyright (C) 2003, Petr Krysl
                 
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/
#ifndef floatarray_h
#   define floatarray_h

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

namespace OOOFS {

class IntArray ; 
class FloatMatrix ;

/**
	Class representing vector of real numbers. This array can grow or shrink to
	desired dimension. The lower value index of array is 1, 
	upper depends on array size. 
 */
class FloatArray
{

public: 

  FloatArray (int = 0) ;                                 
  /**
  Copy constructor. Creates the array from another array.
  */
  FloatArray (const FloatArray&);                         

  virtual ~FloatArray () { if(values) delete [] values; } 

  FloatArray&  operator=  (const FloatArray&);  

#ifdef DEBUG
  /** Coefficient access function. Returns value of coeffiicient at given 
  position of the receiver. Provides 1-based indexing access.
  @param i position of coefficient in array
  */
  double&      at (int i) ;
  /** Coefficient access function. Returns l-value of coeffiicient at given 
  position of the receiver. Provides 1-based indexing access.
  @param i position of coefficient in array
  */
  double       at (int i) const;
#else
  inline double&     at (int i) { return values[i-1] ;}
  inline double      at (int i) const { return values[i-1] ;}
#endif

  /**
  Coefficient access function. Returns value of coeffiicient at given 
  position of the receiver. Provides 0-based indexing access.
  @param i position of coefficient in array
  */
  double&      operator()(int i)
  {
#       ifdef DEBUG
    assert(i < size);
#       endif
    return values[i];
  }
  /**
  Coefficient access function. Returns value of coeffiicient at given 
  position of the receiver. Provides 0-based indexing access.
  @param i position of coefficient in array
  */
  const  double&       operator()(int i) const 
  {
#       ifdef DEBUG
    assert(i < size);
#       endif
    return values[i];
  }

  /** Checks size of receiver towards requested bounds.
  Current implementation will call exit(1), if dimension
  mismatch found.
  @param i required size of receiver
  */
  void         checkBounds (int) const;
  /**
  Checks size of receiver towards values stored in loc array.
  (Expands the receiver if loc points to coefficients beyond the size of receiver).
  */
  void         checkSizeTowards (const IntArray& loc) ;
  /** Checks size of receiver towards requested bounds.
  If dimension mismatch, size is adjusted accordingly.
  The new coefficients are initialized to zero;
  @param allocChunk if reallocation needed, an aditional space for allocChunk values 
  */
  void         resize (int, int allocChunk = 0) ;
  /** Resizes the size of the receiver to requested bounds. Memory allocation always happens, more preferably use
  resize() function instead. 
  */
  void         hardResize (int);
  /**
  Returns nonzero if all coefficients of the receiver are 0, else returns zero.
  */
  int          containsOnlyZeroes () const;
  /// Returns the size of receiver.
  int          giveSize () const         { return size ;}
  /// Returns nonzero if receiver is not empty.
  int          isNotEmpty () const       { return (size != 0) ;}
  /// Returns nonzero if receiver is empty.
  int          isEmpty () const       { return (size == 0) ;}
  /**
  Switches the sign of every coefficient of receiver. 
  @return receiver.
  */
  FloatArray*  negated () ;
  /**
  Print receiver on stdout. Usefull for debugging.
  */
  void         printYourself () const;
  ///	Zeroes all coeficients of receiver.
  void         zero ();
  /* new reference like membre functions */
  /**
  Receiver becomes the result of the product of aMatrix and anArray.
  Adjusts the size of receiver if necessary.
  */
  void  beProductOf (const FloatMatrix& aMatrix, const FloatArray& anArray) ;
  /**
  Receiver becomes the result of the product of aMatrix^T and anArray.
  Adjusts the size of receiver if necessary.
  */
  void  beTProductOf (const FloatMatrix& aMatrix, const FloatArray& anArray) ;
  /**
  Adds array src to receiver. If the receiver's size is zero, it adjusts its size
  to size of src array. If recever's size is nonzero and different from src 
  array size an error is generated.
  */
  void  add         (const FloatArray& src);
  /**
  Substracts array src to receiver. If the receiver's size is zero, it adjusts its size
  to size of src array. If recever's size is nonzero and different from src 
  array size an error is generated.
  */
  void  substract   (const FloatArray& src);
  /**
  Extract sub vector form src array and stores the result into receiver.
  @param src source vector for sub vector
  @param indx Determines sub vector. Receiver size will be indx max value,
  and on i-th position of subVector will be src(indx->at(i)) value.
  */
  void  beSubArrayOf (const FloatArray& src, const IntArray& indx) ;
  /**
  Adds the given vector as sub-vector to receiver. The sub-vector values will be added to receivers
  corresponding receiver's values at positions (si,...,si+src.size). The size of receiver will be 
  adjusted, if necesary.
  @param src the sub-vector to be added
  @param si determines the position (receiver's 1-based index) of first src value to be added.
  */
  void addSubVector (const FloatArray& src, int si);
  /**
  Assembles the array fe (typically, the load vector of a finite
  element) into the receiver, using loc as location array.
  @param fe array to be assebled. 
  @param loc array of code numbers) - src(i) value will 
  be added to receiver value at position loc(i) 
  (if this loc(i) value is nonzero).
  */
  void         assemble (const FloatArray& fe, const IntArray& loc) ;
  /**
  Computes vector product of vectors given as parameters (v1 x v2) and stores the
  result into receiver.
  */
  void  beVectorProductOf (const FloatArray& v1, const FloatArray& v2);
  /**
  Computes the distance between position represented by receiver and position given as parameter. 
  */
  double       distance (const FloatArray& ) const;
  /**
  Returns the receiver 'a' rotated according the change-of-base matrix r.
  @param r Rotation matrix.
  @param mode If mode == 't' the method performs the operation  a = t(transp) * r,
  else if mode = 'n' then the method performs the operation  a = t * r .
  @return modified receiver.
  */
  void  rotatedWith (FloatMatrix& r, char mode) ;
  /* old pointer like member functions */
  FloatArray*  add (FloatArray*) ;
  FloatArray*  substract (FloatArray*);
  FloatArray*  Add (FloatArray* b) {return this->GiveCopy()->add(b);}
  FloatArray*  Substract (FloatArray* b) {return this->GiveCopy()->substract(b);}
  FloatArray*  GiveCopy () const              { return this->Times(1.) ;}
  FloatArray*  GiveSubArray (IntArray*) ;
  double*      givePointer ()  const     { return values ;} // see above
  FloatArray*  beCopyOf (FloatArray*) ;
  FloatArray*  setValuesToZero ();
  FloatArray*  rotatedWith (FloatMatrix* r, char mode) ;
  FloatArray*  times (double) ;
  FloatArray*  Times (double) const ;
  double       distance (FloatArray* );
  FloatArray*  VectorProduct (FloatArray*);
  /**
  Normalizes receiver. Eucleidian norm is used, after operation receiver 
  will have this norm equal to 1.0.
  @return modified receiver
  */
  FloatArray*  normalize ();

  friend double  dotProduct (double*,double*,int) ;
  /**
  Returns the dot product of the first i coefficienst of the two 
  arrays p1 and p2. 
  @return value of the dot product
  */
  friend double  dotProduct (const FloatArray& p1, const FloatArray& p2, int i);

protected:

  int      size;
  int allocatedSize ;
  double*  values ;

} ;

}

#endif








