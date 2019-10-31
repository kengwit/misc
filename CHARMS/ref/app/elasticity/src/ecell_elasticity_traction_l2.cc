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
#include "ecell_elasticity_traction_l2.h"
#include "evalpt.h"
#include "square_fixed_matrix.h"
#include "gcell_line_l2.h"

const char ECELL_ELASTICITY_TRACTION_L2::_type_name[] = "elasticity_traction_l2";

ECELL_ELASTICITY_TRACTION_L2::ECELL_ELASTICITY_TRACTION_L2 (GCELL *gcell)
  : ECELL_ELASTICITY::ECELL_ELASTICITY (gcell)
{
}

bool ECELL_ELASTICITY_TRACTION_L2::assemble_stiffness_matrix (FIELD_VECTOR *u,
                                                     LES *les, FIELD_VECTOR *geometry)
{
  // No stiffness
  return true;
}

bool ECELL_ELASTICITY_TRACTION_L2::assemble_tractions_load (FIELD_VECTOR *u,
                                                      LES *les, FIELD_VECTOR *geometry)
{
  LOAD_ON_GCELL_TRACTION *load = this->traction_load();
  if (load) {
    GAUSS_1d_integration_rule_t *gr = gauss_rule_1d (2);
    vector <DOFMAPPER_VECTOR *> dv = dofmappers ();
    const int nbfuns = dv.size ();
    for (int qp = 0; qp < gr->npoints; qp++) {
      POINT param_loc (gr->ip[qp].xi, 0, 0);
      EVALPT evalpt (u, gcell (), param_loc);
      evalpt.eval ();
      CHECK (evalpt.nbfuns () == nbfuns, EXCEPTION_BAD_VALUE,;);
      POINT at = geometry->evaluate (geometry, this->gcell(), param_loc);
      FIXED_VECTOR<3> traction = load->traction(at);
      FIXED_VECTOR<3> tangent_vector(0);
      for (int I = 0; I < nbfuns; I++) {
        double Nx_I = evalpt.N_x(I);
        FEN *fen =  evalpt.fen(I);
        BFUN_DOFPARAM_PAIR_ID dpid = geometry->dofparam_id (fen);
        FIELD_PAIR<3> *fp = geometry->field_pair (dpid);
        FIXED_VECTOR<3> v = fp->dofparam ();
        tangent_vector.add (Nx_I, v);
      }
      double fact = gr->ip[qp].Wgt * tangent_vector.l2_norm();
      // now loop over all nodes 
      for (int I = 0; I < nbfuns; I++) {
        double N_I = evalpt.N (I);
        double values[3];
        values[0] = traction(0)*fact*N_I;
        values[1] = traction(1)*fact*N_I;
        values[2] = traction(2)*fact*N_I;
        les->rhs_add (dv[I], values);
      }
    }
    gauss_free_1d (gr);
  }
  return true;
}

bool ECELL_ELASTICITY_TRACTION_L2::assemble_body_load (FIELD_VECTOR *u,
                                                      LES *les, FIELD_VECTOR *geometry)
{
  LOAD_ON_GCELL_BODY *body_load = this->body_load();
  if (body_load) {
    GAUSS_3d_integration_rule_t *gr = gauss_rule_3d (2);
    vector <DOFMAPPER_VECTOR *> dv = dofmappers ();
    const int nbfuns = dv.size ();
    for (int qp = 0; qp < gr->npoints; qp++) {
      POINT param_loc (gr->ip[qp].xi, gr->ip[qp].eta, gr->ip[qp].theta);
      EVALPT evalpt (u, gcell (), param_loc);
      evalpt.eval (geometry);
      CHECK (evalpt.nbfuns () == nbfuns, EXCEPTION_BAD_VALUE,;);
      POINT at = geometry->evaluate (geometry, this->gcell(), param_loc);
      FIXED_VECTOR<3> force_density = body_load->force_density(at);
      double fact = gr->ip[qp].Wgt * evalpt.detJ ();
      // now loop over all nodes 
      for (int I = 0; I < nbfuns; I++) {
        double N_I = evalpt.N (I);
        double values[3];
        values[0] = force_density(0)*fact*N_I;
        values[1] = force_density(1)*fact*N_I;
        values[2] = force_density(2)*fact*N_I;
        les->rhs_add (dv[I], values);
      }
    }
    gauss_free_3d (gr);
  }
  return true;
}

bool ECELL_ELASTICITY_TRACTION_L2::assemble_prescribed_u_load (FIELD_VECTOR *u,
                                                      LES *les, FIELD_VECTOR *geometry)
{
  return true;
}

static ECELL_ELASTICITY *
make (GCELL *gcell)
{
  return (new ECELL_ELASTICITY_TRACTION_L2 (gcell));
}

bool
ECELL_ELASTICITY_TRACTION_L2::register_make_func ()
{
  return ECELL_ELASTICITY::register_make_func (make, string (GCELL_LINE_L2::TYPE_NAME),
                                               string ("traction_l2"));
}

