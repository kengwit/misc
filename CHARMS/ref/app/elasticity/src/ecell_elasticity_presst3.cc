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
#include "ecell_elasticity_presst3.h"
#include "evalpt.h"
#include "square_fixed_matrix.h"
#include "gcell_solid_presst3.h"

ECELL_ELASTICITY_PRESST3::ECELL_ELASTICITY_PRESST3 (GCELL *gcell)
  : ECELL_ELASTICITY::ECELL_ELASTICITY (gcell)
{
}

bool ECELL_ELASTICITY_PRESST3::assemble_stiffness_matrix (FIELD_VECTOR *u,
                                                          LES *les, FIELD_VECTOR *geometry)
{
  return true;
}

bool ECELL_ELASTICITY_PRESST3::assemble_load_terms (FIELD_VECTOR *u,
                                               LES *les, FIELD_VECTOR *geometry)
{
  vector <DOFMAPPER_VECTOR *> dv = dofmappers ();
  const int nbfuns = dv.size ();
  // Loop over all quadrature points
  { // for each integration point
    // Handle Dirichlet boundary conditions
    for (int I = 0; I < nbfuns; I++) {
      bool constrained[3];
      dv[I]->constrained (constrained);
      if (constrained[0] || constrained[1] || constrained[2]) { // need to include the contribution to the rhs
        EVALPT evalpt (u, geometry, gcell (), 0.25, 0.25, 0.25);
        evalpt.eval ();
        CHECK (evalpt.nbfuns () == nbfuns, EXCEPTION_BAD_VALUE,;);
        POINT at = geometry->evaluate (geometry, this->gcell(), 0.25, 0.25, 0.25);
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
        double N_zI = evalpt.N_z (I);
        BFUN_DOFPARAM_PAIR_ID dpid = dv[I]->dofparam_id ();
        FIXED_VECTOR<3> prescribed_u = u->field_pair (dpid)->dofparam ();
        if (! prescribed_u.all_zeros ()) {
          double fact = (-1) * (1. / 6) * evalpt.detJ ();
          for (int J = 0; J < nbfuns; J++) {
            double N_xJ = evalpt.N_x (J);
            double N_yJ = evalpt.N_y (J);
            double N_zJ = evalpt.N_z (J);
            double values[3];
            double K[3][3];
            if (is_isotropic) { // RAW is this order correct (B'*C*B)?
              btcb_3d_iso_symm_C (N_xJ, N_yJ, N_zJ, N_xI, N_yI, N_zI, lambda, mu, fact, K);
            } else {
              btcb_3d_aniso_symm_C (N_xJ, N_yJ, N_zJ, N_xI, N_yI, N_zI, C, fact, K);
            }
            for (int m = 0; m < 3; m++) {
              values[m] = K[m][0] * prescribed_u(0) + K[m][1] * prescribed_u(1) + K[m][2] * prescribed_u(2);
            }
            cerr << "Assembling reaction at " << evalpt.bfun(I)->fen()->id() << " to " << evalpt.bfun(J)->fen()->id() << " " << values[0] << " " << values[1] << " " << values[2] << endl;
            les->rhs_add (dv[J], values);
          }
        } // only if non-zero prescribed displacements
      }
    }
  }
  return true;
}

static ECELL_ELASTICITY *
make (GCELL *gcell)
{
  return (new ECELL_ELASTICITY_PRESST3 (gcell));
}

bool
ECELL_ELASTICITY_PRESST3::register_make_func ()
{
  return ECELL_ELASTICITY::register_make_func (make, string (GCELL_SOLID_T3::TYPE_NAME), string ("default"));
}

