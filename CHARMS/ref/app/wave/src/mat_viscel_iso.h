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
#ifndef MAT_VISCEL_ISO_H
# define MAT_VISCEL_ISO_H

#include "mat_viscel.h"
#include "point.h"
#include "func.h"

class MAT_VISCEL_ISO : public MAT_VISCEL {
  
 public:

  static const string TYPE_NAME;

 public: // class functions /////////////////////////////////////////
  
  static bool register_make_func ();
  
 public: // object methods //////////////////////////////////////////

  MAT_VISCEL_ISO (DB *db, string name);

  string type () const { return string (MAT_VISCEL::TYPE_NAME); }

  /**
     Get material stiffness matrix: versions for 1-, 2-, and 3-D.
  */
  bool mat_stiffness (POINT &at, double C[1][1]);
  bool mat_stiffness (POINT &at, double C[3][3]);
  bool mat_stiffness (POINT &at, double C[6][6]);

  bool is_isotropic () const { return true; }

  double lambda (POINT &at) {
    return _lambda (at)(0);
  }

  double mu (POINT &at) {
    return _mu (at)(0);
  }

  double  rho (POINT &at) {
    return _rho (at)(0);
  }

  double  alpha (POINT &at) {
    return _alpha (at)(0);
  }

  double  viscosity (POINT &at) {
    return _viscosity (at)(0);
  }

  double var (string name, POINT &at) {
    if      (name == "lambda") return lambda (at);
    else if (name == "mu")  return mu  (at);
    else if (name == "rho") return rho (at);
    else if (name == "viscosity") return viscosity (at);
    else {
      EXCEPTION_BAD_ACCESS e; throw e;
    }
  }

 private: // object data ////////////////////////////////////////////

  FUNC<1>                        _lambda;
  FUNC<1>                        _mu;
  FUNC<1>                        _rho;
  FUNC<1>                        _alpha;
  FUNC<1>                        _viscosity;

};

#endif

