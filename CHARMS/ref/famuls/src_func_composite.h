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
#ifndef SRC_FUNC_COMPOSITE_H
# define SRC_FUNC_COMPOSITE_H

#include "src_func.h"
#include "src_func_point.h"

class SRC_FUNC_COMPOSITE: public SRC_FUNC {

 public:

  SRC_FUNC_COMPOSITE (list <SRC_FUNC *> member_funcs) : SRC_FUNC() {
    _member_funcs.clear();
    for (list <SRC_FUNC *>::const_iterator i = member_funcs.begin(); i != member_funcs.end(); i++) {
      _member_funcs.push_back (*i);
    }
  }
  
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
