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
#ifndef MAT_ELASTICITY_ISO_H
# define MAT_ELASTICITY_ISO_H

#include "mat_elasticity.h"
#include "point.h"
#include "func.h"

class MAT_ELASTICITY_ISO : public MAT_ELASTICITY {
  
 public:

  static const string TYPE_NAME;

 public: // class functions /////////////////////////////////////////
  
  static bool register_make_func ();
  
 public: // object methods //////////////////////////////////////////

  MAT_ELASTICITY_ISO (DB *db, string name);

  string type () const { return string (MAT_ELASTICITY::TYPE_NAME); }

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

  double sigma11 (POINT &at) {
    return _sigma11 (at)(0);
  }
  
   double sigma22 (POINT &at) {
    return _sigma22 (at)(0);
  }

  double sigma33  (POINT &at) {
    return _sigma33 (at)(0);
  }

  double sigma12 (POINT &at) {
    return _sigma12 (at)(0);
  }
  
  double  sigma13 (POINT &at) {
    return _sigma13 (at)(0);
  }
  double  sigma32 (POINT &at) {
    return _sigma32 (at)(0);
  }
  
  
  
  
  double var (string name, POINT &at) {
    if      (name == "lambda") return lambda (at);
    else if (name == "mu")  return mu  (at);
    else if (name == "rho") return rho (at);
    else if (name == "alpha") return rho (at);
    else if (name == "sigma11") return sigma11 (at);
    else if (name == "sigma22") return sigma22 (at);
    else if (name == "sigma33") return sigma33 (at);
    else if (name == "sigma12") return sigma12 (at);
    else if (name == "sigma13") return sigma13 (at);
    else if (name == "sigma32") return sigma32 (at);
    else {
      EXCEPTION_BAD_ACCESS e; throw e;
    }
  }

 private: // object data ////////////////////////////////////////////

  FUNC<1>                        _lambda;
  FUNC<1>                        _mu;
  FUNC<1>                        _rho;
  FUNC<1>                        _alpha;
  FUNC<1>                        _sigma11;
  FUNC<1>                        _sigma22;
  FUNC<1>                        _sigma33;
  FUNC<1>                        _sigma12;
  FUNC<1>                        _sigma13;
  FUNC<1>                        _sigma32;
};

#endif

