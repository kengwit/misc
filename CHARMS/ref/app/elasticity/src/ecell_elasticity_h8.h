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
#ifndef ECELL_ELASTICITY_H8_H
# define ECELL_ELASTICITY_H8_H

#include "ecell.h"
#include "ecell_elasticity.h"
#include "field.h"
#include "proto_elasticity.h"
#include "gaussP.h"
#include "evalpt.h"
#include "func.h"

class ECELL_ELASTICITY_H8 : public ECELL_ELASTICITY {

 public: // class functions /////////////////////////////////////////

  static bool register_make_func ();

 public: // object functions ////////////////////////////////////////

  ECELL_ELASTICITY_H8 (GCELL *gcell);

  bool assemble_stiffness_matrix (FIELD_VECTOR *u, LES *les, FIELD_VECTOR *geometry);

  bool assemble_prescribed_u_load (FIELD_VECTOR *u, LES *les, FIELD_VECTOR *geometry);
  bool assemble_body_load (FIELD_VECTOR *u, LES *les, FIELD_VECTOR *geometry);
  bool assemble_tractions_load (FIELD_VECTOR *u, LES *les, FIELD_VECTOR *geometry);

 private: // class data /////////////////////////////////////////////

  static const char _type_name[];
  
};

#endif
