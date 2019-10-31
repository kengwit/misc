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
#include "ecell_steady_diffusion_t10.h"
#include "evalpt.h"
#include "square_fixed_matrix.h"
#include "gcell_solid_t10.h"

const char ECELL_STEADY_DIFFUSION_T10::TYPE_NAME[] = "steady_diffusion_t10";

typedef struct gp_coord_t {
  double r, s, t, u;
}              gp_coord_t;

static const gp_coord_t gps[4] = {
  { 0.13819660,0.13819660,0.13819660,0.58541020 },
  { 0.58541020,0.13819660,0.13819660,0.13819660 },
  { 0.13819660,0.58541020,0.13819660,0.13819660 },
  { 0.13819660,0.13819660,0.58541020,0.13819660 },
};
/* The following factor includes W_l=1/4 and c=1/6;
   see Hughes, The finite element method, 1987. */
static const double GPW = 0.041666666666666666667 /* = (1.0/4/6) */;

ECELL_STEADY_DIFFUSION_T10::ECELL_STEADY_DIFFUSION_T10 (GCELL *gcell)
  : ECELL_STEADY_DIFFUSION::ECELL_STEADY_DIFFUSION (gcell)
{
}

bool ECELL_STEADY_DIFFUSION_T10::assemble_conductivity_matrix (FIELD_SCALAR *phi,
                                                              LES *les, FIELD_VECTOR *geometry)
{
  vector <DOFMAPPER_SCALAR *> dv = dofmappers ();
  const int nbfuns = dv.size ();
  SQUARE_FIXED_MATRIX k (nbfuns);
  // First, compute the conductivity matrix
  k.zero ();
  for (int qp = 0; qp < 4; qp++) { // for each integration point
    POINT param_loc (gps[qp].r, gps[qp].s, gps[qp].t);
    EVALPT evalpt (phi, gcell (), param_loc);
    evalpt.eval (geometry);
    CHECK (evalpt.nbfuns () == nbfuns, EXCEPTION_BAD_VALUE,;);
    POINT at = geometry->evaluate (geometry, this->gcell(), param_loc);
    double fact = mat()->conductivity(at) * GPW * evalpt.detJ ();
    // now loop over all nodes 
    for (int I = 0; I < nbfuns; I++) {
      double N_xI = evalpt.N_x (I);
      double N_yI = evalpt.N_y (I);
      double N_zI = evalpt.N_z (I);
      for (int J = 0; J < nbfuns; J++) {
        double N_xJ = evalpt.N_x (J);
        double N_yJ = evalpt.N_y (J);
        double N_zJ = evalpt.N_z (J);
        k(I,J) += (N_xI * N_xJ + N_yI * N_yJ + N_zI * N_zJ) * fact;
      }
    }
  }
  // Finally, we're ready to assemble
  for (unsigned int I = 0; I < k.n(); I++) {
    for (unsigned int J = 0; J < k.n(); J++) {
      double values[1][1] = { { k(I,J) } };
      les->lhs_add (dv[I], dv[J], values);
    }
  }
  return true;
}

bool ECELL_STEADY_DIFFUSION_T10::assemble_source_terms (FIELD_SCALAR *phi,
                                                       LES *les, FIELD_VECTOR *geometry)
{
  vector <DOFMAPPER_SCALAR *> dv = dofmappers ();
  const int nbfuns = dv.size ();
  // Loop over all quadrature points
  for (int qp = 0; qp < 4; qp++) { // for each integration point
    POINT param_loc (gps[qp].r, gps[qp].s, gps[qp].t);
    EVALPT evalpt (phi, gcell (), param_loc);
    evalpt.eval (geometry);
    CHECK (evalpt.nbfuns () == nbfuns, EXCEPTION_BAD_VALUE,;);
    POINT at = geometry->evaluate (geometry, this->gcell(), param_loc);
    // Handle Dirichlet boundary conditions
    for (int I = 0; I < nbfuns; I++) {
      bool constrained[1];
      dv[I]->constrained (constrained);
      if (constrained[0]) {       // need to include the contribution to the rhs
        double N_xI = evalpt.N_x (I);
        double N_yI = evalpt.N_y (I);
        double N_zI = evalpt.N_z (I);
        BFUN_DOFPARAM_PAIR_ID dpid = dv[I]->dofparam_id ();
        FIXED_VECTOR<1> prescribed_phi = phi->field_pair (dpid)->dofparam ();
        double fact = (-1) * mat()->conductivity(at) * prescribed_phi[0] * GPW * evalpt.detJ ();
        for (int J = 0; J < nbfuns; J++) {
          double N_xJ = evalpt.N_x (J);
          double N_yJ = evalpt.N_y (J);
          double N_zJ = evalpt.N_z (J);
          double values[1];
          values[0] = (N_xI * N_xJ + N_yI * N_yJ  + N_zI * N_zJ) * fact;
          les->rhs_add (dv[J], values);
        }
      }
    }
    // Handle internal heat generation
    double fact = mat()->internal_heat_generation_density(at) * GPW * evalpt.detJ ();
    if (fact != 0) {
      for (int I = 0; I < nbfuns; I++) {
        double N = evalpt.N (I);
        double values[1];
        values[0] = N * fact;
        les->rhs_add (dv[I], values);
      }
    }
  }
  return true;
}

static ECELL_STEADY_DIFFUSION *
make (GCELL *gcell)
{
  return (new ECELL_STEADY_DIFFUSION_T10 (gcell));
}

bool
ECELL_STEADY_DIFFUSION_T10::register_make_func ()
{
  return ECELL_STEADY_DIFFUSION::register_make_func (make,
                                                     string (GCELL_SOLID_T10::TYPE_NAME),
                                                     string ("default"));
}

