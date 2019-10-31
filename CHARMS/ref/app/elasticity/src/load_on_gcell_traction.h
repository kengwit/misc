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
#ifndef LOAD_ON_GCELL_TRACTION_H
# define LOAD_ON_GCELL_TRACTION_H

#include "load_on_gcell.h"
#include "point.h"
#include "func.h"

class LOAD_ON_GCELL_TRACTION : public LOAD_ON_GCELL {
  
 public:

  static const string TYPE_NAME;

 public: // class functions /////////////////////////////////////////
  
  static bool register_make_func ();
  
 public: // object methods //////////////////////////////////////////

  LOAD_ON_GCELL_TRACTION (DB *db, string name);

  string type () const { return string (LOAD_ON_GCELL::TYPE_NAME); }

  POINT traction (POINT &at) {
    return _traction (at);
  }
 
  POINT traction (POINT4 &at) {
    return _traction (at);
  }
 
  FUNC<3> &func() { return _traction;}

  double var (string name, POINT &at) {
    if      (name == "traction_x") return traction (at)(0);
    else if (name == "traction_y") return traction (at)(1);
    else if (name == "traction_z") return traction (at)(2);
    else {
      EXCEPTION_BAD_ACCESS e; throw e;
    }
  }

  double var (string name, POINT4 &at) {
    if      (name == "traction_x") return traction (at)(0);
    else if (name == "traction_y") return traction (at)(1);
    else if (name == "traction_z") return traction (at)(2);
    else {
      EXCEPTION_BAD_ACCESS e; throw e;
    }
  }

 private: // object data ////////////////////////////////////////////

  FUNC<3>                        _traction;

};

#endif

