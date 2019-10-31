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
#ifndef ECELL_MODL_H20_H
# define ECELL_MODL_H20_H

#include "ecell.h"
#include "ecell_modl.h"
#include "field.h"
#include "evs_ooofs.h"
#include "proto_modl.h"
#include "gaussP.h"
#include "evalpt.h"
#include "func.h"

class ECELL_MODL_H20 : public ECELL_MODL {

 public: // class functions /////////////////////////////////////////

  static bool register_make_func ();

 public: // object functions ////////////////////////////////////////

  ECELL_MODL_H20 (GCELL *gcell);

  bool assemble_stiffness_matrix (FIELD_VECTOR *u, EVS_OOOFS *evs, FIELD_VECTOR *geometry);
  bool assem_consistent_mass( FIELD_VECTOR *u, EVS_OOOFS *evs, FIELD_VECTOR *geometry, double mmult);
  bool assemble_geom_stiffness_matrix (FIELD_VECTOR *u, EVS_OOOFS *evs, FIELD_VECTOR *geometry,double p);
  bool assemble_geom_stiffness_matrix (FIELD_VECTOR *u, EVS_OOOFS *evs, FIELD_VECTOR *geometry);
  double strain_energy (FIELD_VECTOR *u, FIELD_VECTOR *geometry);
  bool calc_stress (POINT &param_loc,FIELD_VECTOR *u, FIELD_VECTOR *geometry ,double stress[6], double strain[6],double lambda, double mu);
  double mass (FIELD_VECTOR *geometry);

 private: // class data /////////////////////////////////////////////

  static const char _type_name[];
  
};

#endif
