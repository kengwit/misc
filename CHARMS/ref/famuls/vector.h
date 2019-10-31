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
#ifndef VECTOR_H
#   define VECTOR_H


#include <cstdio>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include "famexception.h"

template <typename FLOAT_TYPE>
class VECTOR {
  
 public: // class declarations //////////////////////////////////////
  
  typedef FLOAT_TYPE elem_t;
  typedef elem_t* ptr_t;
  typedef const elem_t* const_ptr_t;
  typedef elem_t* iterator_t;
  typedef const elem_t* const_iterator_t;
  typedef elem_t& ref_t;
  typedef const elem_t& const_ref_t;
  
  /**
     Produces same-size vector, with all components copied.
  */
  VECTOR (const VECTOR<FLOAT_TYPE>& V);
  
  /**
     Size n, all elements set to zero.
   */
  VECTOR (const unsigned int n);
  
  virtual ~VECTOR ();
  
  /**
     Set all elements to s.
   */
  VECTOR<FLOAT_TYPE> & operator= (const FLOAT_TYPE s);
  
  /**
     Copy all elements.
   */
  VECTOR<FLOAT_TYPE> & operator= (const VECTOR<FLOAT_TYPE>& V);
  
  /**
     Copy all elements for different-type vectors.
   */
  template<typename FLOAT_TYPE2>
    VECTOR<FLOAT_TYPE> & operator= (const VECTOR<FLOAT_TYPE2>& V);
  
  /**
     Scalar product.
   */
  FLOAT_TYPE operator* (const VECTOR<FLOAT_TYPE>& V) const;
  
  /**
   * Return square of the $l_2$-norm.
   */
  FLOAT_TYPE l2_norm_squared () const;
  
  /**
   Return the mean value of the elements.
   */
  FLOAT_TYPE mean_value () const;
  
  /**
   Return the $l_1$-norm of the vector, i.e. the sum of the absolute values.
   */
  FLOAT_TYPE l1_norm () const;
  
  /**
   Return the $l_2$-norm of the vector, i.e. the square root of the sum of the
   squares of the elements.
   */
  FLOAT_TYPE l2_norm () const;
  
  /**
   Return the maximum absolute value of the
   elements of this vector, which is the
   L_infinity-norm of a vector.
   */
  FLOAT_TYPE linfinity_norm () const;
  
  
  /**
     Resize the vector to N.  If dont_reinitialize==true, the vector
     is not initialized to zero; otherwise it is.
   */ 
  void resize (const unsigned int N,
               const bool dont_reinitialize=false);

  /**
     Set to zero.
  */
  void zero ();
  /**
     Size of vector.
   */
  unsigned int size () const { return _size; }
  
  iterator_t begin ();
  
  const_iterator_t begin () const;
  
  iterator_t end ();
  
  const_iterator_t end () const;  
  
  /**
     Access element.
  */
  FLOAT_TYPE operator() (const unsigned int i) const;
  
  /**
   Access element as writable.
   */
  FLOAT_TYPE& operator() (const unsigned int i);
  
  /**
     Add U += V.
  */
  VECTOR<FLOAT_TYPE> & operator += (const VECTOR<FLOAT_TYPE> &V);
  
  /**
     Subtract U -= V.
   */
  VECTOR<FLOAT_TYPE> & operator -= (const VECTOR<FLOAT_TYPE> &V);
  /**
     Add a scalar.
   */
  void add (const FLOAT_TYPE s);
  
  /**
    VECTOR addition.
   */
  void add (const VECTOR<FLOAT_TYPE>& V);
  
  /**
   Addition U+=a*V.
   */
  void add (const FLOAT_TYPE a, const VECTOR<FLOAT_TYPE>& V);
  
  /**
   Addition U+=a*V+b*W.
   */
  void add (const FLOAT_TYPE a, const VECTOR<FLOAT_TYPE>& V,
            const FLOAT_TYPE b, const VECTOR<FLOAT_TYPE>& W);
  
  /**
   Addition U=s*U+V.
   */
  void sadd (const FLOAT_TYPE s, const VECTOR<FLOAT_TYPE>& V);
  
  /**
   Addition U=s*U+a*V.
   */
  void sadd (const FLOAT_TYPE s, const FLOAT_TYPE a, const VECTOR<FLOAT_TYPE>& V);
  
  /**
   Addition U=s*U+a*V+b*W.
   */
  void sadd (const FLOAT_TYPE s, const FLOAT_TYPE a,
             const VECTOR<FLOAT_TYPE>& V, const FLOAT_TYPE b, const VECTOR<FLOAT_TYPE>& W);
  
  /**
   Addition U=s*U+a*V+b*W+c*X.
   */
  void sadd (const FLOAT_TYPE s, const FLOAT_TYPE a,
             const VECTOR<FLOAT_TYPE>& V, const FLOAT_TYPE b, const VECTOR<FLOAT_TYPE>& W, 
             const FLOAT_TYPE c, const VECTOR<FLOAT_TYPE>& X);
  
  /**
   Scale each element.
   */
  void scale (const FLOAT_TYPE factor);
  
  /**
   Assignment  U=a*V. 
   */
  void assign (const FLOAT_TYPE a, const VECTOR<FLOAT_TYPE>& V);
  
  /**
   Assignment U=a*V+b*W.
   */
  void assign (const FLOAT_TYPE a, const VECTOR<FLOAT_TYPE>& V,
            const FLOAT_TYPE b, const VECTOR<FLOAT_TYPE>& W);
  
 protected: // object data //////////////////////////////////////////
  
  unsigned int  _size;
  unsigned int  _dim;
  FLOAT_TYPE   *_data;

};

template <typename FLOAT_TYPE>
inline
typename VECTOR<FLOAT_TYPE>::iterator_t VECTOR<FLOAT_TYPE>::begin () {
  return &_data[0];
}


template <typename FLOAT_TYPE>
inline
typename VECTOR<FLOAT_TYPE>::const_iterator_t VECTOR<FLOAT_TYPE>::begin () const {
  return &_data[0];
}


template <typename FLOAT_TYPE>
inline
typename VECTOR<FLOAT_TYPE>::iterator_t VECTOR<FLOAT_TYPE>::end () {
  return &_data[_size];
}


template <typename FLOAT_TYPE>
inline
typename VECTOR<FLOAT_TYPE>::const_iterator_t VECTOR<FLOAT_TYPE>::end () const {
  return &_data[_size];
}


template <typename FLOAT_TYPE>
static inline FLOAT_TYPE squared (const FLOAT_TYPE x) {
  return x*x;
}


template <typename FLOAT_TYPE>
VECTOR<FLOAT_TYPE>::VECTOR (const unsigned int n) : _size(0), _dim(0), _data(0)
{
  resize (n, false);
}


template <typename FLOAT_TYPE>
VECTOR<FLOAT_TYPE>::VECTOR (const VECTOR<FLOAT_TYPE>& v) : _size(v.size()), _dim(v.size()), _data(0)
{
  if (_size) {
    _data = new FLOAT_TYPE[_dim];
    CHECK (_data != 0, EXCEPTION_NULL_PTR,;);
    copy (v.begin(), v.end(), begin());
  }
}


template <typename FLOAT_TYPE>
void VECTOR<FLOAT_TYPE>::resize (const unsigned int n, const bool dont_reinitialize) {
  if (n==0) {
    if (_data) delete[] _data;
    _data = 0;
    _dim = _size = 0;
    return;
  };
  
  if (n>_dim) {
    if (_data) delete[] _data;
    _data = new FLOAT_TYPE[n];
    CHECK (_data != 0, EXCEPTION_NULL_PTR,;);
    _dim = n;
  };
  _size = n;
  if (!dont_reinitialize)
    zero ();
}


template <typename FLOAT_TYPE>
VECTOR<FLOAT_TYPE>::~VECTOR ()
{
  if (_data) delete[] _data;
}


template <typename FLOAT_TYPE>
void VECTOR<FLOAT_TYPE>::zero () {
  if (_size>0)
    fill (begin(), end(), 0.);
}

template <typename FLOAT_TYPE>
FLOAT_TYPE VECTOR<FLOAT_TYPE>::operator * (const VECTOR<FLOAT_TYPE>& v) const
{
  if (&v == this)
    return l2_norm_squared();
  
  CHECK (_size == v._size, EXCEPTION_BAD_ACCESS,;);
  
  FLOAT_TYPE sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;
  
  const_iterator_t ptr  = begin(),
    vptr = v.begin(),
    eptr = ptr + (_size/4)*4;
  while (ptr!=eptr) {
    sum0 += (*ptr++ * *vptr++);
    sum1 += (*ptr++ * *vptr++);
    sum2 += (*ptr++ * *vptr++);
    sum3 += (*ptr++ * *vptr++);
  };
  while (ptr != end())
    sum0 += *ptr++ * *vptr++;
  
  return sum0+sum1+sum2+sum3;
}


template <typename FLOAT_TYPE>
FLOAT_TYPE VECTOR<FLOAT_TYPE>::l2_norm_squared () const
{
  FLOAT_TYPE sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;
  
  const_iterator_t ptr  = begin(),
    eptr = ptr + (_size/4)*4;
  while (ptr!=eptr) {
    sum0 += squared(*ptr++);
    sum1 += squared(*ptr++);
    sum2 += squared(*ptr++);
    sum3 += squared(*ptr++);
  };
  while (ptr != end())
    sum0 += squared(*ptr++);
  
  return sum0+sum1+sum2+sum3;
}


template <typename FLOAT_TYPE>
FLOAT_TYPE VECTOR<FLOAT_TYPE>::mean_value () const
{
  FLOAT_TYPE sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0; 
  
  const_iterator_t ptr  = begin(),
    eptr = ptr + (_size/4)*4;
  while (ptr!=eptr) {
    sum0 += *ptr++;
    sum1 += *ptr++;
    sum2 += *ptr++;
    sum3 += *ptr++;
  };
  while (ptr != end())
    sum0 += *ptr++;
  
  return (sum0+sum1+sum2+sum3)/size();
}


template <typename FLOAT_TYPE>
FLOAT_TYPE VECTOR<FLOAT_TYPE>::l1_norm () const
{
  FLOAT_TYPE sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;

  const_iterator_t ptr  = begin(),
    eptr = ptr + (_size/4)*4;
  while (ptr!=eptr) {
    sum0 += fabs(*ptr++);
    sum1 += fabs(*ptr++);
    sum2 += fabs(*ptr++);
    sum3 += fabs(*ptr++);
  };
  while (ptr != end())
    sum0 += fabs(*ptr++);
  
  return sum0+sum1+sum2+sum3;
}


template <typename FLOAT_TYPE>
FLOAT_TYPE VECTOR<FLOAT_TYPE>::l2_norm () const
{
  return sqrt(l2_norm_squared());
}


template <typename FLOAT_TYPE>
FLOAT_TYPE VECTOR<FLOAT_TYPE>::linfinity_norm () const {
  FLOAT_TYPE max0=0., max1=0., max2=0., max3=0.;
  for (unsigned int i=0; i<(_size/4); ++i) {
    if (max0<fabs(_data[4*i]))   max0=fabs(_data[4*i]);
    if (max1<fabs(_data[4*i+1])) max1=fabs(_data[4*i+1]);
    if (max2<fabs(_data[4*i+2])) max2=fabs(_data[4*i+2]);
    if (max3<fabs(_data[4*i+3])) max3=fabs(_data[4*i+3]);
  };
  for (unsigned int i=(_size/4)*4; i<_size; ++i)
    if (max0<fabs(_data[i]))
      max0 = fabs(_data[i]);
  
  return max (max(max0, max1),
	      max(max2, max3));
}


template <typename FLOAT_TYPE>
VECTOR<FLOAT_TYPE>& VECTOR<FLOAT_TYPE>::operator += (const VECTOR<FLOAT_TYPE>& v)
{
  add (v);
  return *this;
}


template <typename FLOAT_TYPE>
VECTOR<FLOAT_TYPE>& VECTOR<FLOAT_TYPE>::operator -= (const VECTOR<FLOAT_TYPE>& v)
{
  CHECK (_size == v._size, EXCEPTION_BAD_ACCESS,;);
  
  iterator_t i_ptr = begin(),
    i_end = end();
  const_iterator_t v_ptr = v.begin();
  while (i_ptr!=i_end)
    *i_ptr++ -= *v_ptr++;
  
  return *this;
}


template <typename FLOAT_TYPE>
void VECTOR<FLOAT_TYPE>::add (const FLOAT_TYPE v)
{
  iterator_t i_ptr = begin(),
    i_end = end();
  while (i_ptr!=i_end)
    *i_ptr++ += v;
}


template <typename FLOAT_TYPE>
void VECTOR<FLOAT_TYPE>::add (const VECTOR<FLOAT_TYPE>& v)
{
  CHECK (_size == v._size, EXCEPTION_BAD_ACCESS,;);
  
  iterator_t i_ptr = begin(),
    i_end = end();
  const_iterator_t v_ptr = v.begin();
  while (i_ptr!=i_end)
    *i_ptr++ += *v_ptr++;
}


template <typename FLOAT_TYPE>
void VECTOR<FLOAT_TYPE>::add (const FLOAT_TYPE a, const VECTOR<FLOAT_TYPE>& v)
{
  CHECK (_size == v._size, EXCEPTION_BAD_ACCESS,;);
  
  iterator_t i_ptr = begin(),
    i_end = end();
  const_iterator_t v_ptr = v.begin();
  while (i_ptr!=i_end)
    *i_ptr++ += a * *v_ptr++;
}


template <typename FLOAT_TYPE>
void VECTOR<FLOAT_TYPE>::add (const FLOAT_TYPE a, const VECTOR<FLOAT_TYPE>& v,
                              const FLOAT_TYPE b, const VECTOR<FLOAT_TYPE>& w)
{
  CHECK (_size == v._size, EXCEPTION_BAD_ACCESS,;);
  CHECK (_size == w._size, EXCEPTION_BAD_ACCESS,;);
  iterator_t i_ptr = begin(),
    i_end = end();
  const_iterator_t v_ptr = v.begin(),
    w_ptr = w.begin();
  while (i_ptr!=i_end)
    *i_ptr++ += a * *v_ptr++ + b * *w_ptr++;
}


template <typename FLOAT_TYPE>
void VECTOR<FLOAT_TYPE>::sadd (const FLOAT_TYPE x, const VECTOR<FLOAT_TYPE>& v)
{
  CHECK (_size == v._size, EXCEPTION_BAD_ACCESS,;);
  iterator_t i_ptr = begin(),
    i_end = end();
  const_iterator_t v_ptr = v.begin();
  for (; i_ptr!=i_end; ++i_ptr)
    *i_ptr = x * *i_ptr  + *v_ptr++;
}


template <typename FLOAT_TYPE>
void VECTOR<FLOAT_TYPE>::sadd (const FLOAT_TYPE x, const FLOAT_TYPE a, const VECTOR<FLOAT_TYPE>& v)
{
  CHECK (_size == v._size, EXCEPTION_BAD_ACCESS,;);
  iterator_t i_ptr = begin(),
    i_end = end();
  const_iterator_t v_ptr = v.begin();
  for (; i_ptr!=i_end; ++i_ptr)
    *i_ptr = x * *i_ptr  +  a * *v_ptr++;
}


template <typename FLOAT_TYPE>
void VECTOR<FLOAT_TYPE>::sadd (const FLOAT_TYPE x, const FLOAT_TYPE a,
                               const VECTOR<FLOAT_TYPE>& v,
                               const FLOAT_TYPE b, const VECTOR<FLOAT_TYPE>& w)
{
  CHECK (_size == v._size, EXCEPTION_BAD_ACCESS,;);
  CHECK (_size == w._size, EXCEPTION_BAD_ACCESS,;);
  iterator_t i_ptr = begin(),
    i_end = end();
  const_iterator_t v_ptr = v.begin(),
    w_ptr = w.begin();
  for (; i_ptr!=i_end; ++i_ptr)
    *i_ptr = x * *i_ptr  +  a * *v_ptr++  + b * *w_ptr++;
}


template <typename FLOAT_TYPE>
void VECTOR<FLOAT_TYPE>::sadd (const FLOAT_TYPE x, const FLOAT_TYPE a,
                               const VECTOR<FLOAT_TYPE>& v, const FLOAT_TYPE b,
                               const VECTOR<FLOAT_TYPE>& w, const FLOAT_TYPE c, const VECTOR<FLOAT_TYPE>& y)
{
  CHECK (_size == v._size, EXCEPTION_BAD_ACCESS,;);
  CHECK (_size == w._size, EXCEPTION_BAD_ACCESS,;);
  CHECK (_size == y._size, EXCEPTION_BAD_ACCESS,;);

  iterator_t i_ptr = begin(),
    i_end = end();
  const_iterator_t v_ptr = v.begin(),
    w_ptr = w.begin(),
    y_ptr = y.begin();
  
  for (; i_ptr!=i_end; ++i_ptr)
    *i_ptr = (x * *i_ptr)  +  (a * *v_ptr++)  +  (b * *w_ptr++)  + (c * *y_ptr++);
}


template <typename FLOAT_TYPE>
void VECTOR<FLOAT_TYPE>::scale (const FLOAT_TYPE factor)
{
  iterator_t ptr=begin(), eptr=end();
  while (ptr!=eptr)
    *ptr++ *= factor;
}


template <typename FLOAT_TYPE>
void VECTOR<FLOAT_TYPE>::assign (const FLOAT_TYPE a, const VECTOR<FLOAT_TYPE>& u,
                              const FLOAT_TYPE b, const VECTOR<FLOAT_TYPE>& v)
{
  CHECK (_size == u._size, EXCEPTION_BAD_ACCESS,;);
  CHECK (_size == v._size, EXCEPTION_BAD_ACCESS,;);
  
  iterator_t i_ptr = begin(),
    i_end = end();
  const_iterator_t u_ptr = u.begin(),
    v_ptr = v.begin();
  while (i_ptr!=i_end)
    *i_ptr++ = a * *u_ptr++  + b * *v_ptr++;
}


template <typename FLOAT_TYPE>
void VECTOR<FLOAT_TYPE>::assign (const FLOAT_TYPE a, const VECTOR<FLOAT_TYPE>& u)
{
  CHECK (_size == u._size, EXCEPTION_BAD_ACCESS,;);
  iterator_t i_ptr = begin(),
    i_end = end();
  const_iterator_t u_ptr = u.begin();
  while (i_ptr!=i_end)
    *i_ptr++ = a * *u_ptr++;
}

template <typename FLOAT_TYPE>
VECTOR<FLOAT_TYPE>& VECTOR<FLOAT_TYPE>::operator = (const FLOAT_TYPE s)
{
  fill (begin(), end(), s);
  return *this;
}


template <typename FLOAT_TYPE>
VECTOR<FLOAT_TYPE>&
VECTOR<FLOAT_TYPE>::operator = (const VECTOR<FLOAT_TYPE>& v)
{
  if (v._size != _size)
    resize (v._size, true);
  if (_size!=0)
    copy (v.begin(), v.end(), begin());
  
  return *this;
}


template <typename FLOAT_TYPE>
template<typename FLOAT_TYPE2>
VECTOR<FLOAT_TYPE>&
VECTOR<FLOAT_TYPE>::operator = (const VECTOR<FLOAT_TYPE2>& v)
{
  if (v.size() != _size)
    resize (v.size(), true);
  if (_size!=0)
    copy (v.begin(), v.end(), begin());
  
  return *this;
}


#endif
