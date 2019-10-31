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
#ifndef MAT_DIFFUSION_ANISO_H
# define MAT_DIFFUSION_ANISO_H

#include "mat.h"
#include "point.h"
#include "square_fixed_matrix.h"
#include "func.h"
#include "mat_diffusion.h"

class MAT_DIFFUSION_ANISO : public MAT_DIFFUSION {
  
 public: // class data //////////////////////////////////////////////

  static const string TYPE_NAME;

 public: // class functions /////////////////////////////////////////
  
  static bool register_make_func ();
  
 public: // object methods //////////////////////////////////////////

  MAT_DIFFUSION_ANISO (DB *db, string name);

  string type_name () const { return string (TYPE_NAME); }

  double conductivity (POINT &at, int dir) {
    return _conductivity[dir](at)(0);
  }
  double conductivity_r (POINT &at) {
    return conductivity (at, 0);
  }
  double conductivity_s (POINT &at) {
    return conductivity (at, 1);
  }
  double conductivity_t (POINT &at) {
    return conductivity (at, 2);
  }
  double conductivity (POINT &at) {
    return conductivity (at, 0);
  }
  void conductivity (POINT &at, SQUARE_FIXED_MATRIX &kappa) {
    kappa.zero();
    if (kappa.n() == 3) {
      const int ndims = 3;
      for (int dir = 0; dir < ndims; dir++) {
        double k = conductivity (at, dir);
        FIXED_VECTOR<3> d = _dir[dir](at);
        for (int i = 0; i < ndims; i++) {
          for (int j = 0; j < ndims; j++) {
            kappa(i,j) += k * d(i)*d(j);
          }
        }
      }
    } else if (kappa.n() == 2) {
      const int ndims = 2;
      for (int dir = 0; dir < ndims; dir++) {
        double k = conductivity (at, dir);
        FIXED_VECTOR<3> d = _dir[dir](at);
        for (int i = 0; i < ndims; i++) {
          for (int j = 0; j < ndims; j++) {
            kappa(i,j) += k * d(i)*d(j);
          }
        }
      }
    } else if (kappa.n() == 1) {
      kappa(0,0) = _conductivity[0](at)(0);
    }
  }

  double internal_heat_generation_density (POINT &at) {
    return _internal_heat_generation_density (at)(0);
  }

  double var (string name, POINT &at) {
    if      (name == "conductivity_r") return conductivity_r (at);
    else if (name == "internal_heat_generation_density") return internal_heat_generation_density (at);
    else {
      EXCEPTION_BAD_ACCESS e; throw e;
    }
  }

 private: // object data ////////////////////////////////////////////

  FUNC<3> _dir[3]; // direction in which the conductivity matrix is diagonal
  FUNC<1> _conductivity[3];
  FUNC<1> _internal_heat_generation_density;

};

#endif

