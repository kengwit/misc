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
#include "ecell_errest_steady_diffusion_l2.h"
#include "gcell_line_l2.h"
#include "evalpt.h"
#include "proto_errest.h"

const char ECELL_ERREST_STEADY_DIFFUSION_L2::TYPE_NAME[] = "errest_steady_diffusion_l2";

ECELL_ERREST_STEADY_DIFFUSION_L2::ECELL_ERREST_STEADY_DIFFUSION_L2 (GCELL *gcell)
  : ECELL_ERREST<1>::ECELL_ERREST (gcell) 
{
}

void
ECELL_ERREST_STEADY_DIFFUSION_L2::assemble_error (FIELD_SCALAR *field,
                                 FIELD_VECTOR *geometry, PROTO_ERREST<1> *proto)
{
  POINT param_loc (0, 0, 0);
  EVALPT evalpt (field, gcell (), param_loc);
  evalpt.eval (geometry);
  const size_t nbfuns = evalpt.nbfuns (); 
  FIXED_VECTOR<1> grad;
  grad = 0; double area = evalpt.detJ();
  for (size_t i = 0; i < nbfuns; i++) {
    BFUN *bfun = evalpt.bfun (i);
    FEN *fen = bfun->fen ();
    BFUN_DOFPARAM_PAIR_ID dpid = field->dofparam_id (fen);
    FIELD_PAIR<1> *fp = field->field_pair (dpid);
    FIXED_VECTOR<1> dofparam = fp->dofparam ();
    double N_x = evalpt.N_x (i);
    double fact = N_x * evalpt.detJ();
    grad.add(fact, dofparam);
  }
  // Error ~ magnitude of the gradient
  double gradnorm = grad.l2_norm ();
  // Now assemble the error 
  for (size_t i = 0; i < nbfuns; i++) {
    BFUN *bfun = evalpt.bfun (i);
    proto->accum (bfun, gradnorm, area);
  }  
}

static ECELL_ERREST<1> *
make (GCELL *gcell)
{
  return (new ECELL_ERREST_STEADY_DIFFUSION_L2 (gcell));
}

bool
ECELL_ERREST_STEADY_DIFFUSION_L2::register_make_func ()
{
  return ECELL_ERREST<1>::register_make_func (make, string (GCELL_LINE_L2::TYPE_NAME), string ("default"));
}
