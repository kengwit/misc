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
#include "ecell_errest_elasticity_plane_strain_q4.h"
#include "gcell_surf_q4.h"
#include "evalpt.h"
#include "proto_errest.h"
#include "mat_map.h"

const char ECELL_ERREST_ELASTICITY_PLANE_STRAIN_Q4::TYPE_NAME[] = "errest_elasticity_plane_strain_q4";

ECELL_ERREST_ELASTICITY_PLANE_STRAIN_Q4::ECELL_ERREST_ELASTICITY_PLANE_STRAIN_Q4 (GCELL *gcell)
  : ECELL_ERREST<3>::ECELL_ERREST (gcell) 
{
  _mat = 0;
}

#include "btcb.h"

void
ECELL_ERREST_ELASTICITY_PLANE_STRAIN_Q4::assemble_error (FIELD_VECTOR *u,
                                 FIELD_VECTOR *geometry, class PROTO_ERREST<3> *proto)
{
  // Compute the stiffness matrix submatrix linking two nodes
  { // one-point quadrature
    double err = 0;
    POINT param_loc (0.0, 0.0, 0.0);
    EVALPT evalpt (u, gcell (), param_loc);
    evalpt.eval (geometry);
    POINT at = geometry->evaluate (geometry, this->gcell(), param_loc);
    double C[6][6];
    double lambda = 0;
    double mu = 0;
    bool is_isotropic = _mat->is_isotropic();
    if (is_isotropic) {
      lambda = _mat->var("lambda", at);
      mu = _mat->var ("mu", at);
    } else {
      _mat->mat_stiffness (at, C);
    }
    double vol = 4 * evalpt.detJ (); // one-point quadrature
    // now loop over all nodes 
    size_t nbfuns = evalpt.nbfuns();
    for (size_t I = 0; I < nbfuns; I++) {
      double N_xI = evalpt.N_x (I);
      double N_yI = evalpt.N_y (I);
      double N_zI = 0;
      BFUN *bfunI = evalpt.bfun (I);
      FEN *fenI = bfunI->fen ();
      BFUN_DOFPARAM_PAIR_ID dpidI = u->dofparam_id (fenI);
      FIELD_PAIR<3> *fpI = u->field_pair (dpidI);
      FIXED_VECTOR<3> uI = fpI->dofparam ();
      for (size_t J = I; J < nbfuns; J++) {
        double N_xJ = evalpt.N_x (J);
        double N_yJ = evalpt.N_y (J);
        double N_zJ = 0;
        BFUN *bfunJ = evalpt.bfun (J);
        FEN *fenJ = bfunJ->fen ();
        BFUN_DOFPARAM_PAIR_ID dpidJ = u->dofparam_id (fenJ);
        FIELD_PAIR<3> *fpJ = u->field_pair (dpidJ);
        FIXED_VECTOR<3> uJ = fpJ->dofparam ();
        double K[3][3];
        if (is_isotropic) {
          btcb_3d_iso_symm_C (N_xJ, N_yJ, N_zJ, N_xI, N_yI, N_zI, lambda, mu, vol, K);
        } else {
          btcb_3d_aniso_symm_C (N_xJ, N_yJ, N_zJ, N_xI, N_yI, N_zI, C, vol, K);
        }
        double derr = (uJ(0)*(K[0][0]*uI(0)+
                              K[0][1]*uI(1)+
                              K[0][2]*uI(2))+
                       uJ(1)*(K[1][0]*uI(0)+
                              K[1][1]*uI(1)+
                              K[1][2]*uI(2))+
                       uJ(2)*(K[2][0]*uI(0)+
                              K[2][1]*uI(1)+
                              K[2][2]*uI(2)));
        err += derr;
        if (I != J) { // also for the lower triangle
          err += derr;
        }
      }
    }
    // Error ~ magnitude of the gradient
    for (size_t I = 0; I < nbfuns; I++) {
      BFUN *bfun = evalpt.bfun (I);
      proto->accum (bfun, err, vol);
    }
  }
}

static ECELL_ERREST<3> *
make (GCELL *gcell)
{
  return (new ECELL_ERREST_ELASTICITY_PLANE_STRAIN_Q4 (gcell));
}

bool
ECELL_ERREST_ELASTICITY_PLANE_STRAIN_Q4::register_make_func ()
{
  return ECELL_ERREST<3>::register_make_func (make, string (GCELL_SURF_Q4::TYPE_NAME), string ("default"));
}


void
ECELL_ERREST_ELASTICITY_PLANE_STRAIN_Q4::attach (class PROTO_ERREST<3> *proto, string ggpath)
{
  ALGO *algo = proto->algo();
  DB *db = algo->mgr()->db();
  // Material
  string matname = db->DB_GET_STRING (ggpath + "/material");
  string mattype = db->DB_GET_STRING ("materials/" + matname + "/type");
  MAT_MGR *mat_mgr = algo->mgr()->mat_mgr();
  MAT_MAP <MAT_ELASTICITY> mat_map (mat_mgr);
  _mat = mat_map.mat (mattype, matname);
  //cerr << "done in ECELL_ERREST_ELASTICITY_PLANE_STRAIN_Q4::attach" << endl;
}
