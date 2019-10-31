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
#ifndef ECELL_ELASTICITY_PRESST3_H
# define ECELL_ELASTICITY_PRESST3_H

#include "ecell.h"
#include "ecell_elasticity.h"
#include "field.h"
#include "les.h"
#include "proto_elasticity.h"
#include "evalpt.h"
#include "func.h"

class ECELL_ELASTICITY_PRESST3 : public ECELL_ELASTICITY {

 public: // class functions /////////////////////////////////////////

  static bool register_make_func ();

 public: // object functions ////////////////////////////////////////

  ECELL_ELASTICITY_PRESST3 (GCELL *gcell);

  bool assemble_stiffness_matrix (FIELD_VECTOR *u, LES *les, FIELD_VECTOR *geometry);

  bool assemble_load_terms (FIELD_VECTOR *u, LES *les, FIELD_VECTOR *geometry);

 private: // class data /////////////////////////////////////////////

  static const char _type_name[] = "elasticity_presst3";
  
};

#endif
