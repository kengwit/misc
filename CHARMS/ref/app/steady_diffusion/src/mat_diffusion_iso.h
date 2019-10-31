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
#ifndef MAT_DIFFUSION_ISO_H
# define MAT_DIFFUSION_ISO_H

#include "mat.h"
#include "point.h"
#include "func.h"
#include "mat_diffusion.h"

class MAT_DIFFUSION_ISO : public MAT_DIFFUSION {
  
 public: // class data //////////////////////////////////////////////

  static const string TYPE_NAME;

 public: // class functions /////////////////////////////////////////
  
  static bool register_make_func ();
  
 public: // object methods //////////////////////////////////////////

  MAT_DIFFUSION_ISO (DB *db, string name);

  string type_name () const { return string (TYPE_NAME); }

  double conductivity (POINT &at) {
    return _conductivity (at)(0);
  }

  double internal_heat_generation_density (POINT &at) {
    return _internal_heat_generation_density (at)(0);
  }

  double var (string name, POINT &at) {
    if      (name == "conductivity") return conductivity (at);
    else if (name == "internal_heat_generation_density") return internal_heat_generation_density (at);
    else {
      EXCEPTION_BAD_ACCESS e; throw e;
    }
  }

 private: // object data ////////////////////////////////////////////

  FUNC<1>                        _conductivity;
  FUNC<1>                        _internal_heat_generation_density;

};

#endif

