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
#ifndef intarray_h
#   define intarray_h

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

namespace OOOFS {

/**
Class implementing an array of integers. Array can grow or shrink to desired dimension.
The lower value index of array is 1, upper depends on array size.
*/
class IntArray
{
  /*
  This class implements an array of integers.
  DESCRIPTION :
  An IntArray stores its coefficients in an array 'value' of size 'size'.
  TASKS :
  - storing and reterning coefficients (method 'at') ;
  - appending another IntArray to itself ;
  - introduced allocatedSize variable to allow dynamic rescaling of array
  size possibly without memory realocation. At startup array occupies space
  given by allocatedSpace = size. Then there can be 
  1) further request for resizeing array to smaller dimension
  then we only change size wariable, but allocatedSize
  variable remain untouched - expecting possible array grow and then re-using
  previously allocated space.
  2) if further request for growing then is necessary memory realocation.
  This process is controlled in resize member function. 
  REMARK :
  see Remark 2 in file "floatarry.hxx".
  */




public:
  /// Constructor for zero sized array
  IntArray (int=0) ;                                   // constructor
  /** Copy constructor. Creates the array from another array.
  */
  IntArray (const IntArray&);                        // copy constructor
  /// Destructor.
  ~IntArray ()  { if(values) delete [] values ;}      // destructor

  /// Assingnment operator
  IntArray&  operator=  (const IntArray&);              // assignment: cleanup and copy


#     ifdef DEBUG
  /** Coefficient access function. Returns l-value of coeffiicient at given 
  position of the receiver.
  @param i position of coefficient in array
  */
  int&    at (int i) ; 
  /** Coefficient access function. Returns value of coeffiicient at given 
  position of the receiver.
  @param i position of coefficient in array
  */
  int     at (int i) const;
#     else
  int&    at (int i)                  { return values[i-1] ;}
  int     at (int i) const            { return values[i-1] ;}
#     endif
  /**
  Coefficient access function. Returns value of coeffiicient at given 
  position of the receiver. Provides 0-based indexing access.
  @param i position of coefficient in array
  */
  int&      operator()(int i)
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
  const  int&       operator()(int i) const 
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
  void       checkBounds (int i) const;
  /** Checks size of receiver towards requested bounds.
  If dimension mismatch, size is adjusted accordingly.
  Warning: after this operation array values are in undefined state, programmer should
  zero receiver 
  @param allocChunk if reallocation needed, an aditional space for allocChunk values 
  */
  void       resize (int n, int allocChunk = 0) ;
  /**
  Appends array b at the end of receiver.
  @param b array to be appended at the end of receiver
  @param allocChunk if reallocation needed, an aditional space for allocChunk values 
  will be allocated to prevent excessive realocation
  */
  void       followedBy (const IntArray& b, int allocChunk = 0) ;
  /**
  Appends given Number at the end of receiver.
  @param b value to be appended
  @param allocChunk if reallocation needed, an aditional space for allocChunk values 
  will be allocated to prevent excessive realocation
  */
  void       followedBy (const int b, int allocChunk = 0);
  /// Returns the size of receiver.
  int        giveSize () const           { return size ;}
  /// Checks if receiver is empty (i.e., zero sized).
  int        isEmpty ()  const           { return size==0 ;}
  /**
  Finds index of first occurence of given value in array. If such value is not presented,
  returns zero value.
  @param value scanned value
  @return index of first value in array, otherwise zero
  */
  int        findFirstIndexOf (int value)  const ;
  /// Prints receiver on stdin.
  void       printYourself () const ;
  /// Sets all component to zero.
  void       zero() ;

private:

  /// size of array
  int   size;
  /// allocated size for array
  int  allocatedSize ;
  /// stored values
  int*  values ;

} ;


template <class operation> int  
quickSortPartition (IntArray& arry, int l, int r, operation op) {
  int i=l-1, j=r;
  int v = arry.at(r);
  int swap;

  for (;;) {
    while (( op(arry.at(++i), v)) < 0 );
    while (( op(v, arry.at(--j))) < 0 ) if (j==l) break;
    if (i >= j) break;
    swap = arry.at(i); arry.at(i) = arry.at(j); arry.at(j) = swap;
  }
  swap = arry.at(i); arry.at(i) = arry.at(r); arry.at(r) = swap;
  return i;
}



template <class operation> void quickSort (IntArray& arry, int l, int r, operation op) {
  if (r<=l) return;
  int i = quickSortPartition (arry, l, r, op);
  quickSort (arry, l, i-1, op);
  quickSort (arry, i+1, r, op);
}


/** 
Sorts the receiver using quiksort algorithm.
@param op is Function object, required to have member function int class::operator() (int, int),
must return a negative value if first argument is less than the second,
zero if the arguments are equal, and a positive number otherwise.
*/
template <class operation> void sort (IntArray& arry, operation op)	{quickSort (arry, 1, arry.giveSize(), op);}

}


#endif


