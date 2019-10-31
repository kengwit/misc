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
#include "ecell_modl_h8_bbar.h"
#include "evalpt.h"
#include "square_fixed_matrix.h"
#include "gcell_solid_h8.h"

const char ECELL_MODL_H8_BBAR::_type_name[] = "modl_h8_bbar";

ECELL_MODL_H8_BBAR::ECELL_MODL_H8_BBAR (GCELL *gcell)
  : ECELL_MODL::ECELL_MODL (gcell)
{
}

#include "btcbbar.h"

bool ECELL_MODL_H8_BBAR::assemble_stiffness_matrix (FIELD_VECTOR *u,
                                                     EVS_OOOFS *evs, FIELD_VECTOR *geometry)
{
  GAUSS_3d_integration_rule_t *gr = gauss_rule_3d (2);
  vector <DOFMAPPER_VECTOR *> dv = dofmappers ();
  const int nbfuns = dv.size ();
  double K[3][3];
  // calculate the mean dilatation 
  double N_x_avg=0, N_y_avg=0, N_z_avg=0;
  double elvol=0;
  for (int qp = 0; qp < gr->npoints; qp++) {
    POINT param_loc (gr->ip[qp].xi, gr->ip[qp].eta, gr->ip[qp].theta);
    EVALPT evalpt (u, gcell (), param_loc);
    evalpt.eval (geometry);
    CHECK (evalpt.nbfuns () == nbfuns, EXCEPTION_BAD_VALUE,;);
    double fact = gr->ip[qp].Wgt * evalpt.detJ ();
    elvol += fact;
    for (int I = 0; I < nbfuns; I++) {
      double N_xI = evalpt.N_x (I);
      double N_yI = evalpt.N_y (I);
      double N_zI = evalpt.N_z (I);
      N_x_avg += fact * N_xI;
      N_y_avg += fact * N_yI;
      N_z_avg += fact * N_zI;
    }
  }
  N_x_avg /= elvol;
  N_y_avg /= elvol;
  N_z_avg /= elvol;
  // Compute the stiffness matrix submatrix linking two nodes
  for (int qp = 0; qp < gr->npoints; qp++) {
    POINT param_loc (gr->ip[qp].xi, gr->ip[qp].eta, gr->ip[qp].theta);
    EVALPT evalpt (u, gcell (), param_loc);
    evalpt.eval (geometry);
    CHECK (evalpt.nbfuns () == nbfuns, EXCEPTION_BAD_VALUE,;);
    POINT at = geometry->evaluate (geometry, this->gcell(), param_loc);
    double lambda = 0;
    double mu = 0;
    lambda = mat()->var("lambda", at);
    mu = mat()->var ("mu", at);
    double fact = gr->ip[qp].Wgt * evalpt.detJ ();
    // now loop over all nodes 
    for (int I = 0; I < nbfuns; I++) {
      double N_xI = evalpt.N_x (I);
      double N_yI = evalpt.N_y (I);
      double N_zI = evalpt.N_z (I);
      for (int J = I; J < nbfuns; J++) {
        double N_xJ = evalpt.N_x (J);
        double N_yJ = evalpt.N_y (J);
        double N_zJ = evalpt.N_z (J);
	btcb_3d_iso_symm_C_bbar (N_xJ, N_yJ, N_zJ, N_xI, N_yI, N_zI, N_x_avg, N_y_avg, N_z_avg, lambda, mu, fact, K);
	evs->add_to_op (dv[J], dv[I], K);
        if (I != J) {
          transp (K);
          evs->add_to_op (dv[I], dv[J], K);
        }
      }
    }
  }
  gauss_free_3d (gr);
  return true;
}


bool ECELL_MODL_H8_BBAR::assemble_geom_stiffness_matrix (FIELD_VECTOR *u, EVS_OOOFS *evs, FIELD_VECTOR *geometry,double p)
{
  //This stiffness matrix assumes that material is homogenous isotropic
  GAUSS_3d_integration_rule_t *gr = gauss_rule_3d (2);
  vector <DOFMAPPER_VECTOR *> dv = dofmappers ();
  const int nbfuns = dv.size ();
  double lambda = 0;
  double mu = 0;
  int    qp = 0;
  POINT param_loc (gr->ip[0].xi, gr->ip[0].eta, gr->ip[0].theta);
  POINT at = geometry->evaluate (geometry, this->gcell(), param_loc);
  lambda = mat()->var("lambda", at);
  mu = mat()->var ("mu", at);
  for (qp = 0; qp < 4; qp++) { // for each integration point
    POINT param_loc (gr->ip[qp].xi, gr->ip[qp].eta, gr->ip[qp].theta);
    EVALPT evalpt (u, gcell (), param_loc);
    evalpt.eval (geometry);
    CHECK (evalpt.nbfuns () == nbfuns, EXCEPTION_BAD_VALUE,;);
    POINT at = geometry->evaluate (geometry, this->gcell(), param_loc);
    double fact = gr->ip[qp].Wgt * evalpt.detJ ();
    // now loop over all nodes 
    for (int I = 0; I < nbfuns; I++) {
      double N_xI = evalpt.N_x (I);
      double N_yI = evalpt.N_y (I);
      double N_zI = evalpt.N_z (I);
      for (int J = I; J < nbfuns; J++) {
        double N_xJ = evalpt.N_x (J);
        double N_yJ = evalpt.N_y (J);
        double N_zJ = evalpt.N_z (J);
        double k[3][3] = {{(N_xI*N_xJ+N_yI*N_yJ+N_zI*N_zJ)*p*fact,   0,                                  0},
                           {0,   (N_xI*N_xJ+N_yI*N_yJ+N_zI*N_zJ)*p*fact,                                  0},
                           {0,                                   0,  (N_xI*N_xJ+N_yI*N_yJ+N_zI*N_zJ)*p*fact}};
        evs->add_to_op (dv[J], dv[I], k);
        if (I != J) {
          transp (k);
          evs->add_to_op (dv[I], dv[J], k);
        }
      }
    }
  }
  return true;
}


bool
ECELL_MODL_H8_BBAR::assem_consistent_mass (FIELD_VECTOR *u,EVS_OOOFS *evs, FIELD_VECTOR *geometry, double mmult)
{
  vector <DOFMAPPER_VECTOR *> dv = dofmappers ();
  const int nbfuns = dv.size ();
  GAUSS_3d_integration_rule_t *gr = gauss_rule_3d (2);
  POINT param_loc (gr->ip[0].xi, gr->ip[0].eta, gr->ip[0].theta);
  POINT at = geometry->evaluate (geometry, gcell(), param_loc);
  double rho = mat()->var("rho", at);
  double m[3][3]={{0,0,0},{0,0,0},{0,0,0}};
  for (int qp = 0; qp < gr->npoints; qp++) {
    POINT param_loc (gr->ip[qp].xi, gr->ip[qp].eta, gr->ip[qp].theta);
    EVALPT evalpt (u, gcell(), param_loc);
    evalpt.eval (geometry);
    double fact = rho * gr->ip[qp].Wgt * evalpt.detJ ();
    // now loop over all nodes
    for (int I = 0; I < nbfuns; I++) {
      double N_I = evalpt.N (I);
      for (int J = 0; J < nbfuns; J++) {
        double N_J = evalpt.N (J);
        m[0][0] = m[1][1] = m[2][2] = mmult*fact*N_I*N_J;
        evs->add_to_op (dv[J], dv[I], m);
      }
    }
  }
  return true;
}


static ECELL_MODL *
make (GCELL *gcell)
{
  return (new ECELL_MODL_H8_BBAR (gcell));
}

bool
ECELL_MODL_H8_BBAR::register_make_func ()
{
  return ECELL_MODL::register_make_func (make, string (GCELL_SOLID_H8::TYPE_NAME),
					       string ("bbar"));
}

