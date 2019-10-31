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
#ifndef LOAD_ON_GCELL_BODY_H
# define LOAD_ON_GCELL_BODY_H

#include "load_on_gcell.h"
#include "point.h"
#include "func.h"

class LOAD_ON_GCELL_BODY : public LOAD_ON_GCELL {
  
 public:

  static const string TYPE_NAME;

 public: // class functions /////////////////////////////////////////
  
  static bool register_make_func ();
  
 public: // object methods //////////////////////////////////////////

  LOAD_ON_GCELL_BODY (DB *db, string name);

  string type () const { return string (LOAD_ON_GCELL::TYPE_NAME); }

  POINT force_density (POINT &at) {
    return _force_density (at);
  }

  POINT force_density (POINT4 &at) {
    return _force_density (at);
  }

  double var (string name, POINT &at) {
    if      (name == "force_density_x") return force_density (at)(0);
    else if (name == "force_density_y") return force_density (at)(1);
    else if (name == "force_density_z") return force_density (at)(2);
    else {
      EXCEPTION_BAD_ACCESS e; throw e;
    }
  }

  double var (string name, POINT4 &at) {
    if      (name == "force_density_x") return force_density (at)(0);
    else if (name == "force_density_y") return force_density (at)(1);
    else if (name == "force_density_z") return force_density (at)(2);
    else {
      EXCEPTION_BAD_ACCESS e; throw e;
    }
  }

 private: // object data ////////////////////////////////////////////

  FUNC<3>                        _force_density;

};

#endif

