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
#ifndef FIXED_VECTOR_H
#   define FIXED_VECTOR_H


#include "famuls.h"
#include <cstdio>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include "famexception.h"

template <int DIM>
class FIXED_VECTOR {
  
 public: // class declarations //////////////////////////////////////
  
  /**
     Produces vector with all components initialized to zero.
  */
  FIXED_VECTOR () { this->zero (); }
  /**
     Produces vector with the all components set to the argument.
  */
  FIXED_VECTOR (double s) { *this = s; }
  /**
     Produces vector with all components copied.
  */
  FIXED_VECTOR (const FIXED_VECTOR<DIM>& V);
  /**
     Produces vector with the first two components set to s and t; all others to zero.
  */
  FIXED_VECTOR (const double s, const double t) {
    (*this)(0) = s;
    (*this)(1) = t;
  }
  /**
     Produces vector with the first three components set to s, t, and u; all others to zero.
  */
  FIXED_VECTOR (const double s, const double t, const double u) {
    (*this)(0) = s;
    (*this)(1) = t;
    (*this)(2) = u;
  }
    
  ~FIXED_VECTOR ();
  
  /**
     Get the size of the vector.
  */
  sizet size () const { return DIM; }
  
  /**
     Set all elements to s.
   */
  FIXED_VECTOR<DIM> & operator= (const double s);
  
  /**
     Copy all elements.
   */
  FIXED_VECTOR<DIM> & operator= (const FIXED_VECTOR<DIM>& V);
  
  /**
     Add two vectors.
   */
  FIXED_VECTOR<DIM> operator+ (const FIXED_VECTOR<DIM>& V);
  
  /**
     Subtract one vector from another.
   */
  FIXED_VECTOR<DIM> operator- (const FIXED_VECTOR<DIM>& V);
  
  /**
     Scalar product.
   */
  double operator* (const FIXED_VECTOR<DIM>& V) const;
  /**
    Square of the $l_2$-norm.
   */
  double l2_norm_squared () const;
  
  /**
   Return the mean value of the elements.
   */
  double mean_value () const;
  
  /**
   Return the $l_1$-norm of the fixed_vector, i.e. the sum of the absolute values.
   */
  double l1_norm () const;
  
  /**
   Return the $l_2$-norm of the fixed_vector, i.e. the square root of the sum of the
   squares of the elements.
   */
  double l2_norm () const;
  
  /**
   Return the maximum absolute value of the
   elements of this fixed_vector, which is the L_infinity-norm of a vector.
   */
  double linfinity_norm () const;
  
  /**
     Set to zero.
  */
  void zero ();

  /**
     All zeros?
  */
  bool all_zeros ();
  
  /**
     Access element.
  */
  double operator() (const unsigned int i) const;
  
  /**
   Access element as writable.
   */
  double& operator() (const unsigned int i);
  
  /**
     Access element.
  */
  double operator[] (const unsigned int i) const;
  
  /**
   Access element as writable.
   */
  double& operator[] (const unsigned int i);
  
  /**
     Add U += V.
  */
  FIXED_VECTOR<DIM> & operator += (const FIXED_VECTOR<DIM> &V);
  
  /**
     Subtract U -= V.
   */
  FIXED_VECTOR<DIM> & operator -= (const FIXED_VECTOR<DIM> &V);
  /**
     Add a scalar.
   */
  void add (const double s);
  
  /**
    Vector addition.
   */
  void add (const FIXED_VECTOR<DIM>& V);
  
  /**
   Addition U+=a*V.
   */
  void add (const double a, const FIXED_VECTOR<DIM>& V);
  
  /**
   Addition U+=a*V+b*W.
   */
  void add (const double a, const FIXED_VECTOR<DIM>& V,
            const double b, const FIXED_VECTOR<DIM>& W);
  
  /**
   Addition U=s*U+V.
   */
  void sadd (const double s, const FIXED_VECTOR<DIM>& V);
  
  /**
   Addition U=s*U+a*V.
   */
  void sadd (const double s, const double a, const FIXED_VECTOR<DIM>& V);
  
  /**
   Addition U=s*U+a*V+b*W.
   */
  void sadd (const double s, const double a,
             const FIXED_VECTOR<DIM>& V, const double b, const FIXED_VECTOR<DIM>& W);
  
  /**
   Addition U=s*U+a*V+b*W+c*X.
   */
  void sadd (const double s, const double a,
             const FIXED_VECTOR<DIM>& V, const double b, const FIXED_VECTOR<DIM>& W, 
             const double c, const FIXED_VECTOR<DIM>& X);
  
  /**
   Scale each element.
   */
  FIXED_VECTOR<DIM> & scale (const double factor);
  
  /**
   Assignment  U=a*V. 
   */
  void assign (const double a, const FIXED_VECTOR<DIM>& V);
  
  /**
   Assignment U=a*V+b*W.
   */
  void assign (const double a, const FIXED_VECTOR<DIM>& V,
            const double b, const FIXED_VECTOR<DIM>& W);
  /**
   Assignment U=a*V+b*W+c*Z.
   */
  void assign (const double a, const FIXED_VECTOR<DIM>& V,
            const double b, const FIXED_VECTOR<DIM>& W,
               const double c, const FIXED_VECTOR<DIM>& Z);
  /**
   Assignment U=a*V+b*W+c*Y+d*Z.
   */
  void assign (const double a, const FIXED_VECTOR<DIM>& V,
               const double b, const FIXED_VECTOR<DIM>& W,
               const double c, const FIXED_VECTOR<DIM>& Y,
               const double d, const FIXED_VECTOR<DIM>& Z
               );

  /**
     All components equal?
  */
  bool equal (const FIXED_VECTOR<DIM> &ov);
  
 private: // object data //////////////////////////////////////////
  
  double   _data[DIM];

};


template <int DIM>
inline
double FIXED_VECTOR<DIM>::operator() (const unsigned int i) const
{
  return _data[i];
}


template <int DIM>
inline
double & FIXED_VECTOR<DIM>::operator() (const unsigned int i) 
{
  return _data[i];
}


template <int DIM>
inline
double FIXED_VECTOR<DIM>::operator[] (const unsigned int i) const
{
  return _data[i];
}


template <int DIM>
inline
double & FIXED_VECTOR<DIM>::operator[] (const unsigned int i) 
{
  return _data[i];
}




template <int DIM>
FIXED_VECTOR<DIM>::FIXED_VECTOR (const FIXED_VECTOR<DIM>& v) 
{
  for (int j = 0; j < DIM; j++) _data[j] = v(j);
}

template <int DIM>
FIXED_VECTOR<DIM>::~FIXED_VECTOR () {}


template <int DIM>
void FIXED_VECTOR<DIM>::zero () {
  for (int j = 0; j < DIM; j++) _data[j] = 0;
}


template <int DIM>
bool FIXED_VECTOR<DIM>::all_zeros () {
  for (int j = 0; j < DIM; j++) if (_data[j] != 0) return false;
  return true;
}

template <int DIM>
double FIXED_VECTOR<DIM>::operator * (const FIXED_VECTOR<DIM>& v) const
{
  if (&v == this)
    return l2_norm_squared();
  
  double sum = 0;
  for (int j = 0; j < DIM; j++) sum += _data[j] * v(j);
  return sum;
}

template <int DIM>
double FIXED_VECTOR<DIM>::l2_norm_squared () const
{
  double sum = 0;
  for (int j = 0; j < DIM; j++) sum += (_data[j])*(_data[j]);
  return sum;
}


template <int DIM>
double FIXED_VECTOR<DIM>::mean_value () const
{
  double sum = 0;
  for (int j = 0; j < DIM; j++) sum += _data[j];
  return (sum/DIM);
}

template <int DIM>
double FIXED_VECTOR<DIM>::l1_norm () const
{
  double sum = 0;
  for (int j = 0; j < DIM; j++) sum += fabs(_data[j]);
  return sum;
}


template <int DIM>
double FIXED_VECTOR<DIM>::l2_norm () const
{
  return sqrt(l2_norm_squared());
}


template <int DIM>
double FIXED_VECTOR<DIM>::linfinity_norm () const {
  double max0=0.;
  for (int j = 0; j < DIM; j++) {
    if (max0 < fabs (_data[j])) max0 = fabs (_data[j]);
  }
  return max0;
}

template <int DIM>
FIXED_VECTOR<DIM>& FIXED_VECTOR<DIM>::operator += (const FIXED_VECTOR<DIM>& v)
{
  for (int j = 0; j < DIM; j++) _data[j] += v(j);
  return *this;
}

template <int DIM>
FIXED_VECTOR<DIM> FIXED_VECTOR<DIM>::operator + (const FIXED_VECTOR<DIM>& V)
{
  FIXED_VECTOR<DIM> r = *this;
  return r += V;
}

template <int DIM>
FIXED_VECTOR<DIM> FIXED_VECTOR<DIM>::operator - (const FIXED_VECTOR<DIM>& V)
{
  FIXED_VECTOR<DIM> r = *this;
  return r -= V;
}


template <int DIM>
FIXED_VECTOR<DIM>& FIXED_VECTOR<DIM>::operator -= (const FIXED_VECTOR<DIM>& v)
{
  for (int j = 0; j < DIM; j++) _data[j] -= v(j);
  return *this;
}


template <int DIM>
void FIXED_VECTOR<DIM>::add (const double s)
{
  for (int j = 0; j < DIM; j++) _data[j] += s;
}


template <int DIM>
void FIXED_VECTOR<DIM>::add (const FIXED_VECTOR<DIM>& v)
{
  for (int j = 0; j < DIM; j++) _data[j] += v(j);
}


template <int DIM>
void FIXED_VECTOR<DIM>::add (const double a, const FIXED_VECTOR<DIM>& v)
{
  for (int j = 0; j < DIM; j++) _data[j] += a * v(j);
}


template <int DIM>
void FIXED_VECTOR<DIM>::add (const double a, const FIXED_VECTOR<DIM>& v,
                              const double b, const FIXED_VECTOR<DIM>& w)
{
  for (int j = 0; j < DIM; j++) _data[j] += (a * v(j)) + (b * w(j));
}


template <int DIM>
void FIXED_VECTOR<DIM>::sadd (const double a, const FIXED_VECTOR<DIM>& v)
{
  for (int j = 0; j < DIM; j++) _data[j] = a * _data[j] + v(j);
}


template <int DIM>
void FIXED_VECTOR<DIM>::sadd (const double s, const double a, const FIXED_VECTOR<DIM>& v)
{
  for (int j = 0; j < DIM; j++) _data[j] = s * _data[j] + a * v(j);
}


template <int DIM>
void FIXED_VECTOR<DIM>::sadd (const double s, const double a,
                               const FIXED_VECTOR<DIM>& v,
                               const double b, const FIXED_VECTOR<DIM>& w)
{
  for (int j = 0; j < DIM; j++) _data[j] = s * _data[j] + (a * v(j)) + (b * w(j));
}


template <int DIM>
void FIXED_VECTOR<DIM>::sadd (const double s, const double a,
                               const FIXED_VECTOR<DIM>& v, const double b,
                               const FIXED_VECTOR<DIM>& w, const double c, const FIXED_VECTOR<DIM>& y)
{
   for (int j = 0; j < DIM; j++) _data[j] = s * _data[j] + (a * v(j)) + (b * w(j)) + (c * y(j));
}


template <int DIM>
FIXED_VECTOR<DIM> & FIXED_VECTOR<DIM>::scale (const double factor)
{
  for (int j = 0; j < DIM; j++) _data[j] *= factor;
  return *this;
}


template <int DIM>
void FIXED_VECTOR<DIM>::assign (const double a, const FIXED_VECTOR<DIM>& u,
                                const double b, const FIXED_VECTOR<DIM>& v,
                                const double c, const FIXED_VECTOR<DIM>& z
                                )
{
  for (int j = 0; j < DIM; j++) _data[j] = a * u(j) + b * v(j) + c * z(j);  
}

template <int DIM>
void FIXED_VECTOR<DIM>::assign (const double a, const FIXED_VECTOR<DIM>& u,
                                const double b, const FIXED_VECTOR<DIM>& v,
                                const double c, const FIXED_VECTOR<DIM>& y,
                                const double d, const FIXED_VECTOR<DIM>& z
                                )
{
  for (int j = 0; j < DIM; j++) _data[j] = a * u(j) + b * v(j) + c * y(j) + d * z(j);  
}

template <int DIM>
void FIXED_VECTOR<DIM>::assign (const double a, const FIXED_VECTOR<DIM>& u,
                                const double b, const FIXED_VECTOR<DIM>& v)
{
  for (int j = 0; j < DIM; j++) _data[j] = a * u(j) + b * v(j);  
}


template <int DIM>
void FIXED_VECTOR<DIM>::assign (const double a, const FIXED_VECTOR<DIM>& u)
{
  for (int j = 0; j < DIM; j++) _data[j] = a * u(j);   
}

template <int DIM>
FIXED_VECTOR<DIM>& FIXED_VECTOR<DIM>::operator = (const double s)
{
  for (int j = 0; j < DIM; j++) _data[j] = s;
  return *this;
}


template <int DIM>
FIXED_VECTOR<DIM>&
FIXED_VECTOR<DIM>::operator = (const FIXED_VECTOR<DIM>& v)
{
  for (int j = 0; j < DIM; j++) _data[j] = v(j);  
  return *this;
}

template <int DIM>
bool
FIXED_VECTOR<DIM>::equal (const FIXED_VECTOR<DIM>& ov)
{
  for (int j = 0; j < DIM; j++) if (_data[j] != ov(j)) return false; 
  return true;
}

  /**
     Output.
  */
template <int DIM>
std::ostream & operator << (std::ostream &s, FIXED_VECTOR<DIM> v) {
  s << "(";
  for (sizet j = 0; j < DIM-1; j++) s << v(j) << ",";
  return s << v(DIM-1) << ")";
}

extern FIXED_VECTOR<3> cross_prod( const FIXED_VECTOR<3>& u, const FIXED_VECTOR<3>& v );

#endif


