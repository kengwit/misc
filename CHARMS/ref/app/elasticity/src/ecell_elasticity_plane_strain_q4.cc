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
#include "ecell_elasticity_plane_strain_q4.h"
#include "evalpt.h"
#include "square_fixed_matrix.h"
#include "gcell_surf_q4.h"

const char ECELL_ELASTICITY_PLANE_STRAIN_Q4::_type_name[] = "elasticity_plane_strain_q4";

ECELL_ELASTICITY_PLANE_STRAIN_Q4::ECELL_ELASTICITY_PLANE_STRAIN_Q4 (GCELL *gcell)
  : ECELL_ELASTICITY::ECELL_ELASTICITY (gcell)
{
}

#include "btcb.h"

bool ECELL_ELASTICITY_PLANE_STRAIN_Q4::assemble_stiffness_matrix (FIELD_VECTOR *u,
                                                     LES *les, FIELD_VECTOR *geometry)
{
  GAUSS_2d_integration_rule_t *gr = gauss_rule_2d (2);
  vector <DOFMAPPER_VECTOR *> dv = dofmappers ();
  const int nbfuns = dv.size ();
  double K[3][3];
  // Compute the stiffness matrix submatrix linking two nodes
  for (int qp = 0; qp < gr->npoints; qp++) {
    POINT param_loc (gr->ip[qp].xi, gr->ip[qp].eta, 0);
    EVALPT evalpt (u, gcell (), param_loc);
    evalpt.eval (geometry);
    CHECK (evalpt.nbfuns () == nbfuns, EXCEPTION_BAD_VALUE,;);
    POINT at = geometry->evaluate (geometry, this->gcell(), param_loc);
    double C[6][6];
    double lambda = 0;
    double mu = 0;
    bool is_isotropic = mat()->is_isotropic();
    if (is_isotropic) {
      lambda = mat()->var("lambda", at);
      mu = mat()->var ("mu", at);
    } else {
      mat()->mat_stiffness (at, C);
    }
    double fact = gr->ip[qp].Wgt * evalpt.detJ ();
    // now loop over all nodes 
    for (int I = 0; I < nbfuns; I++) {
      double N_xI = evalpt.N_x (I);
      double N_yI = evalpt.N_y (I);
      double N_zI = 0;
      for (int J = I; J < nbfuns; J++) {
        double N_xJ = evalpt.N_x (J);
        double N_yJ = evalpt.N_y (J);
        double N_zJ = 0;
        if (is_isotropic) {
          btcb_3d_iso_symm_C (N_xJ, N_yJ, N_zJ, N_xI, N_yI, N_zI, lambda, mu, fact, K);
        } else {
          btcb_3d_aniso_symm_C (N_xJ, N_yJ, N_zJ, N_xI, N_yI, N_zI, C, fact, K);
        }
        les->lhs_add (dv[J], dv[I], K);
        if (I != J) {
          transp (K);
          les->lhs_add (dv[I], dv[J], K);
        }
      }
    }
  }
  gauss_free_2d (gr);
  return true;
}

bool ECELL_ELASTICITY_PLANE_STRAIN_Q4::assemble_tractions_load (FIELD_VECTOR *u,
                                                      LES *les, FIELD_VECTOR *geometry)
{
  // No tractions may be applied to this type of element
  return true;
}

bool ECELL_ELASTICITY_PLANE_STRAIN_Q4::assemble_body_load (FIELD_VECTOR *u,
                                                      LES *les, FIELD_VECTOR *geometry)
{
  LOAD_ON_GCELL_BODY *body_load = this->body_load();
  if (body_load) {
    GAUSS_2d_integration_rule_t *gr = gauss_rule_2d (2);
    vector <DOFMAPPER_VECTOR *> dv = dofmappers ();
    const int nbfuns = dv.size ();
    for (int qp = 0; qp < gr->npoints; qp++) {
      POINT param_loc (gr->ip[qp].xi, gr->ip[qp].eta, 0);
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
    gauss_free_2d (gr);
  }
  return true;
}

bool ECELL_ELASTICITY_PLANE_STRAIN_Q4::assemble_prescribed_u_load (FIELD_VECTOR *u,
                                                      LES *les, FIELD_VECTOR *geometry)
{
  GAUSS_2d_integration_rule_t *gr = gauss_rule_2d (2);
  vector <DOFMAPPER_VECTOR *> dv = dofmappers ();
  const int nbfuns = dv.size ();
  // Loop over all quadrature points
  for (int qp = 0; qp < gr->npoints; qp++) {
    // Handle Dirichlet boundary conditions
    for (int I = 0; I < nbfuns; I++) {
      bool constrained[3];
      dv[I]->constrained (constrained);
      if (constrained[0] || constrained[1] || constrained[2]) { // need to include the contribution to the rhs
        POINT param_loc (gr->ip[qp].xi, gr->ip[qp].eta, 0);
        EVALPT evalpt (u, gcell (), param_loc);
        evalpt.eval (geometry);
        CHECK (evalpt.nbfuns () == nbfuns, EXCEPTION_BAD_VALUE,;);
        POINT at = geometry->evaluate (geometry, this->gcell(), param_loc);
        double C[6][6];
        double lambda = 0, mu = 0;
        bool is_isotropic = mat()->is_isotropic();
        if (is_isotropic) {
          lambda = mat()->var("lambda", at);
          mu = mat()->var ("mu", at);
        } else {
          mat()->mat_stiffness (at, C);
        }
        double N_xI = evalpt.N_x (I);
        double N_yI = evalpt.N_y (I);
        double N_zI = 0;
        BFUN_DOFPARAM_PAIR_ID dpid = dv[I]->dofparam_id ();
        FIXED_VECTOR<3> prescribed_u = u->field_pair (dpid)->dofparam ();
        bool constrained[3];
        u->field_pair(dpid)->constrained(constrained);
        if (! prescribed_u.all_zeros ()) {
          double fact = (-1) * gr->ip[qp].Wgt * evalpt.detJ ();
          for (int J = 0; J < nbfuns; J++) {
            double N_xJ = evalpt.N_x (J);
            double N_yJ = evalpt.N_y (J);
            double N_zJ = 0;
            double K[3][3];
            if (is_isotropic) { // RAW is this order correct (B'*C*B)?
              btcb_3d_iso_symm_C (N_xJ, N_yJ, N_zJ, N_xI, N_yI, N_zI, lambda, mu, fact, K);
            } else {
              btcb_3d_aniso_symm_C (N_xJ, N_yJ, N_zJ, N_xI, N_yI, N_zI, C, fact, K);
            }
            double values[3];
            for (int m = 0; m < 3; m++) {
              values[m] = 0;
              for (int n = 0; n < 3; n++) {
                if (constrained[n]) values[m] += K[m][n] * prescribed_u(n);
              }
            }
            les->rhs_add (dv[J], values);
          }
        } // only if non-zero prescribed displacements
      }
    }
  }
  gauss_free_2d (gr);
  return true;
}

static ECELL_ELASTICITY *
make (GCELL *gcell)
{
  return (new ECELL_ELASTICITY_PLANE_STRAIN_Q4 (gcell));
}

bool
ECELL_ELASTICITY_PLANE_STRAIN_Q4::register_make_func ()
{
  return ECELL_ELASTICITY::register_make_func (make, string (GCELL_SURF_Q4::TYPE_NAME), string ("default"));
}

