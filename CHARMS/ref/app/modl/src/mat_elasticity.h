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
#ifndef MAT_ELASTICITY_H
# define MAT_ELASTICITY_H

#include "mat.h"
#include "point.h"
#include "func.h"

class MAT_ELASTICITY : public MAT {
  
 public: // class data //////////////////////////////////////////////

  static const string TYPE_NAME;

 public: // object methods //////////////////////////////////////////

  MAT_ELASTICITY (DB *db, string name) : MAT (name) {}
  virtual ~MAT_ELASTICITY () {}

  string type_name () const { return TYPE_NAME; }

  virtual bool is_isotropic () const { return false; }
  /**
     Get material stiffness matrix: versions for 1-, 2-, and 3-D.
  */
  virtual bool mat_stiffness (POINT &at, double C_mat[1][1]) { return false; }
  virtual bool mat_stiffness (POINT &at, double C_mat[3][3]) { return false; }
  virtual bool mat_stiffness (POINT &at, double C_mat[6][6]) { return false; }

};

#endif

