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
#ifndef SRC_FUNC_POINT_H
# define SRC_FUNC_POINT_H

#include "src_func.h"

class SRC_FUNC_POINT: public SRC_FUNC {

 public:

  SRC_FUNC_POINT (POINT center, double v, double R) : SRC_FUNC(v) { _center = center; _value = v; _R = R; }
  
  ~SRC_FUNC () {}

  FIXED_VECTOR<NUM_COMPONENTS> operator() (POINT& at) {
    POINT r = at - _center;
    double d_I  = r.l2_norm ();
    if (d_I < 1.e-9 * _R) d_I = 1.e-9;
    double a = (_R - d_I) / (_R * d_I);
    FIXED_VECTOR<1> result(_value * a * a);
    return result;
  }

 private:

  double _value;
  POINT  _center;
  double _R;

};

#endif
