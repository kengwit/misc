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
#include "famuls.h"
#include "mat_diffusion_aniso.h"

const string MAT_DIFFUSION_ANISO::TYPE_NAME = "diffusion_aniso";

static MAT *
make (pair <DB *, string> p)
{
  return (new MAT_DIFFUSION_ANISO (p.first, p.second));
}

bool
MAT_DIFFUSION_ANISO::register_make_func ()
{
  return MAT::register_make_func (make, string(TYPE_NAME));
}

MAT_DIFFUSION_ANISO::MAT_DIFFUSION_ANISO (DB *db, string name) : MAT_DIFFUSION (name)
{
  string path = "materials/" + this->name();
  string expr[1];
  expr[0] = db->DB_GET_STRING (path + "/conductivity_r");
  _conductivity[0] = FUNC<1> (expr);
  expr[0] = db->DB_GET_STRING (path + "/dir_r");
  _dir[0] = FUNC<3> (expr);
  //
  expr[0] = db->DB_GET_STRING (path + "/conductivity_s");
  _conductivity[1] = FUNC<1> (expr);
  expr[0] = db->DB_GET_STRING (path + "/dir_s");
  _dir[1] = FUNC<3> (expr);
  //
  expr[0] = db->DB_GET_STRING (path + "/conductivity_t");
  _conductivity[2] = FUNC<1> (expr);
  expr[0] = db->DB_GET_STRING (path + "/dir_t");
  _dir[2] = FUNC<3> (expr);
  //
  expr[0] = db->DB_GET_STRING (path + "/internal_heat_generation_density");
  _internal_heat_generation_density = FUNC<1> (expr);
}
