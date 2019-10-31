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
#ifndef ECELL_wave_L2_H
# define ECELL_wave_L2_H

#include "ecell.h"
#include "ecell_wave.h"
#include "field.h"
#include "les.h"
#include "proto_wave.h"
#include "gaussP.h"
#include "evalpt.h"
#include "func.h"

class ECELL_wave_L2 : public ECELL_wave {

 public: // class functions /////////////////////////////////////////

  static bool register_make_func ();

 public: // object functions ////////////////////////////////////////

  ECELL_wave_L2 (GCELL *gcell);

  bool assemble_stiffness_matrix (FIELD_VECTOR *u, LES *les, FIELD_VECTOR *geometry);

  bool assemble_prescribed_u_load (FIELD_VECTOR *u, LES *les, FIELD_VECTOR *geometry);
  bool assemble_body_load (FIELD_VECTOR *u, LES *les, FIELD_VECTOR *geometry);
  bool assemble_tractions_load (FIELD_VECTOR *u, LES *les, FIELD_VECTOR *geometry);
  FIXED_VECTOR<3>calc_loc ( FIELD_VECTOR *geometry,POINT &param_loc) {};
  vector<double> calc_quants( FIELD_VECTOR *u, FIELD_VECTOR *v,POINT &param_loc,list<string>vars) {};

 private: // class data /////////////////////////////////////////////

  static const char _type_name[];
  
};

#endif
