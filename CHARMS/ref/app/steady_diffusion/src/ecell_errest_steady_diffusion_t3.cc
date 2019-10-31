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
#include "ecell_errest_steady_diffusion_t3.h"
#include "gcell_surf_t3.h"
#include "evalpt.h"
#include "proto_errest.h"

const char ECELL_ERREST_STEADY_DIFFUSION_T3::TYPE_NAME[] = "errest_steady_diffusion_t3";

ECELL_ERREST_STEADY_DIFFUSION_T3::ECELL_ERREST_STEADY_DIFFUSION_T3 (GCELL *gcell)
  : ECELL_ERREST<1>::ECELL_ERREST (gcell) 
{
}

void
ECELL_ERREST_STEADY_DIFFUSION_T3::assemble_error (FIELD_SCALAR *field,
                                 FIELD_VECTOR *geometry, class PROTO_ERREST<1> *proto)
{
  trinumintg_integration_rule_t *tir = trinumintg_integration_rule (1);
  FIXED_VECTOR<1> f_x;
  FIXED_VECTOR<1> f_y;
  f_x = 0; f_y = 0;
  double area = 0;
  for (int qp = 0; qp < tir->npoints; qp++) {
    POINT param_loc (tir->ip[qp].L[0], tir->ip[qp].L[1], 0);
    EVALPT evalpt (field, gcell (), param_loc);
    evalpt.eval (geometry);
    double fact = tir->ip[qp].Wgt * evalpt.detJ ();
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
    }
    // Error ~ magnitude of the gradient
    double gradnorm = sqrt (f_x.l2_norm_squared () + f_y.l2_norm_squared ());
    for (int I = 0; I < nbfuns; I++) {
      BFUN *bfun = evalpt.bfun (I);
      proto->accum (bfun, gradnorm, area);
    }
  }
  
  // and clean up
  trinumintg_free_rule (tir);
}

static ECELL_ERREST<1> *
make (GCELL *gcell)
{
  return (new ECELL_ERREST_STEADY_DIFFUSION_T3 (gcell));
}

bool
ECELL_ERREST_STEADY_DIFFUSION_T3::register_make_func ()
{
  return ECELL_ERREST<1>::register_make_func (make, string (GCELL_SURF_T3::TYPE_NAME), string ("default"));
}
