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
#include "ecell_errest_steady_diffusion_t10.h"
#include "gcell_solid_t10.h"
#include "evalpt.h"
#include "proto_errest.h"

const char ECELL_ERREST_STEADY_DIFFUSION_T10::TYPE_NAME[] = "errest_steady_diffusion_t10";

ECELL_ERREST_STEADY_DIFFUSION_T10::ECELL_ERREST_STEADY_DIFFUSION_T10 (GCELL *gcell)
  : ECELL_ERREST<1>::ECELL_ERREST (gcell) 
{
}

void
ECELL_ERREST_STEADY_DIFFUSION_T10::assemble_error (FIELD_SCALAR *field,
                                 FIELD_VECTOR *geometry, class PROTO_ERREST<1> *proto)
{
  FIXED_VECTOR<1> f_x(0.0), f_y(0.0), f_z(0.0);
  double area = 0;
  { // a single quadrature point
    POINT param_loc (0.25, 0.25, 0.25);
    EVALPT evalpt (field, gcell (), param_loc);
    evalpt.eval (geometry);
    double fact = (1.0 / 6) * evalpt.detJ ();
    area += fact;
    const int nbfuns = evalpt.nbfuns ();
    for (int I = 0; I < nbfuns; I++) {
      BFUN *bfun = evalpt.bfun (I);
      FEN *fen = bfun->fen ();
      BFUN_DOFPARAM_PAIR_ID dpid = field->dofparam_id (fen);
      FIELD_PAIR<1> *fp = field->field_pair (dpid);
      FIXED_VECTOR<1> dofparam = fp->dofparam ();
      double N_x = evalpt.N_x (I);
      f_x.add(fact*N_x, dofparam);
      double N_y = evalpt.N_y (I);
      f_y.add(fact*N_y, dofparam);
      double N_z = evalpt.N_z (I);
      f_z.add(fact*N_z, dofparam);
    }
    // Error ~ magnitude of the gradient
    double gradnorm = sqrt (f_x.l2_norm_squared () + f_y.l2_norm_squared () + f_z.l2_norm_squared ());
    for (int I = 0; I < nbfuns; I++) {
      BFUN *bfun = evalpt.bfun (I);
      proto->accum (bfun, gradnorm, area);
    }
  }
}

static ECELL_ERREST<1> *
make (GCELL *gcell)
{
  return (new ECELL_ERREST_STEADY_DIFFUSION_T10 (gcell));
}

bool
ECELL_ERREST_STEADY_DIFFUSION_T10::register_make_func ()
{
  return ECELL_ERREST<1>::register_make_func (make, string (GCELL_SOLID_T10::TYPE_NAME), string ("default"));
}
