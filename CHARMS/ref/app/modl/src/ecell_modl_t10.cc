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
#include "ecell_modl_t10.h"
#include "evalpt.h"
#include "square_fixed_matrix.h"
#include "gcell_solid_t10.h"

const char ECELL_MODL_T10::TYPE_NAME[] = "modl_t10";

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

/*
  5 point formula
 */
static const gp_coord_t gps5[5] = {
  {0.25, 0.25, 0.25, 0.25 },
  {0.166666666666666666667,0.166666666666666666667,0.166666666666666666667,0.5},
  {0.5,0.166666666666666666667,0.166666666666666666667,0.166666666666666666667},
  {0.166666666666666666667,0.5,0.166666666666666666667,0.166666666666666666667},
  {0.166666666666666666667,0.166666666666666666667,0.5,0.166666666666666666667},            
}; 
static const double GPW5_1 = -0.13333333333333333333; /* -4/5/6  */ 
static const double GPW5_2 =  0.075;  /*9/20/6*/

/*
  16 point formula 
 */

static const gp_coord_t gps16[16] = {
  { .7716429020672371,  .07611903264425430, .07611903264425430 },
  { .07611903264425430, .7716429020672371,  .07611903264425430 },
  { .07611903264425430, .07611903264425430, .7716429020672371  }, 
  { .07611903264425430, .07611903264425430, .07611903264425430 },
  { .1197005277978019,  .07183164526766925, .4042339134672644  },
  { .1197005277978019,  .4042339134672644,  .07183164526766925 }, 
  { .1197005277978019,  .4042339134672644,  .4042339134672644  },
  { .07183164526766925, .1197005277978019,  .4042339134672644  },
  { .07183164526766925, .4042339134672644,  .1197005277978019  },
  { .07183164526766925, .4042339134672644,  .4042339134672644  },
  { .4042339134672644,  .4042339134672644,  .1197005277978019  },
  { .4042339134672644,  .1197005277978019,  .4042339134672644  },
  { .4042339134672644,  .4042339134672644,  .07183164526766925,}, 
  { .4042339134672644,  .07183164526766925, .4042339134672644  },
  { .4042339134672644,  .1197005277978019,  .07183164526766925 },
  { .4042339134672644,  .07183164526766925, .1197005277978019  },
};
static const double GPW16_1 = .05037379410012282;
static const double GPW16_2 = .066542068863329239;  

/*
17 point formula
 */
static const gp_coord_t gps17[17] = {
  {.25, .25, .25, .25 },
  {.7316369079576180, .08945436401,      .08945436401 },
  {.08945436401,      .7316369079576180, .08945436401},
  {.08945436401,      .08945436401,      .7316369079576180},
  {.08945436401,      .08945436401,      .08945436401},
  {.1325810999384657, .02454003792,      .4214394310662522},
  {.1325810999384657, .4214394310662522, .02454003792},
  {.1325810999384657, .4214394310662522, .4214394310662522},
  {.02454003792,      .1325810999384657, .4214394310662522 },
  {.02454003792,      .4214394310662522, .1325810999384657 },
  {.02454003792,      .4214394310662522, .4214394310662522},
  {.4214394310662522, .02454003792,      .4214394310662522},
  {.4214394310662522, .4214394310662522, .02454003792},
  {.4214394310662522, .1325810999384657, .4214394310662522 },
  {.4214394310662522, .4214394310662522, .1325810999384657},
  {.4214394310662522, .1325810999384657, .02454003792 },
  {.4214394310662522, .02454003792,      .1325810999384657},
};

static const double GPW17_1 = .1884185567365411; 
static const double GPW17_2 = .06703858372604275;
static const double GPW17_3 = .04528559236327399;

ECELL_MODL_T10::ECELL_MODL_T10 (GCELL *gcell)
  : ECELL_MODL::ECELL_MODL (gcell)
{
}
#include "btcb.h"

bool ECELL_MODL_T10::assemble_stiffness_matrix (FIELD_VECTOR *u,
                                                              EVS_OOOFS *evs, FIELD_VECTOR *geometry)
{
  //This stiffness matrix assumes that material is homogenous isotropic
  vector <DOFMAPPER_VECTOR *> dv = dofmappers ();
  const int nbfuns = dv.size ();
  // SQUARE_FIXED_MATRIX k (nbfuns);
  double k[3][3];
  // k.zero ();
  double lambda = 0;
  double mu = 0;
  int    qp = 0;
  POINT param_loc (gps[qp].r, gps[qp].s, gps[qp].t);
  POINT at = geometry->evaluate (geometry, this->gcell(), param_loc);
  lambda = mat()->var("lambda", at);
  mu = mat()->var ("mu", at);
  for (qp = 0; qp < 4; qp++) { // for each integration point
    POINT param_loc (gps[qp].r, gps[qp].s, gps[qp].t);
    EVALPT evalpt (u, gcell (), param_loc);
    evalpt.eval (geometry);
    CHECK (evalpt.nbfuns () == nbfuns, EXCEPTION_BAD_VALUE,;);
    POINT at = geometry->evaluate (geometry, this->gcell(), param_loc);
    double fact = GPW* evalpt.detJ ();
    // now loop over all nodes 
    for (int I = 0; I < nbfuns; I++) {
      double N_xI = evalpt.N_x (I);
      double N_yI = evalpt.N_y (I);
      double N_zI = evalpt.N_z (I);
      for (int J = I; J < nbfuns; J++) {
        double N_xJ = evalpt.N_x (J);
        double N_yJ = evalpt.N_y (J);
        double N_zJ = evalpt.N_z (J);
        btcb_3d_iso_symm_C (N_xJ, N_yJ, N_zJ, N_xI, N_yI, N_zI, lambda, mu, fact, k);
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

bool ECELL_MODL_T10::assemble_geom_stiffness_matrix (FIELD_VECTOR *u,
                                                     EVS_OOOFS *evs, FIELD_VECTOR *geometry,
                                                     double p)
{
  //This stiffness matrix assumes that material is homogenous isotropic
  vector <DOFMAPPER_VECTOR *> dv = dofmappers ();
  const int nbfuns = dv.size ();
  double lambda = 0;
  double mu = 0;
  int    qp = 0;
  POINT param_loc (gps[qp].r, gps[qp].s, gps[qp].t);
  POINT at = geometry->evaluate (geometry, this->gcell(), param_loc);
  lambda = mat()->var("lambda", at);
  mu = mat()->var ("mu", at);
  for (qp = 0; qp < 4; qp++) { // for each integration point
    POINT param_loc (gps[qp].r, gps[qp].s, gps[qp].t);
    EVALPT evalpt (u, gcell (), param_loc);
    evalpt.eval (geometry);
    CHECK (evalpt.nbfuns () == nbfuns, EXCEPTION_BAD_VALUE,;);
    POINT at = geometry->evaluate (geometry, this->gcell(), param_loc);
    double fact = GPW* evalpt.detJ ();
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
ECELL_MODL_T10::assem_consistent_mass (FIELD_VECTOR *u,
                                       EVS_OOOFS *evs, FIELD_VECTOR *geometry, double mmult)
{
  vector <DOFMAPPER_VECTOR *> dv = dofmappers ();
  const int nbfuns = dv.size ();
  int qp = 0; double fact; 
  POINT param_loc (gps16[qp].r, gps16[qp].s, gps16[qp].t);
  POINT at = geometry->evaluate (geometry, gcell(), param_loc);
  double rho = mat()->var("rho", at);
  double m[3][3]={{0,0,0},{0,0,0},{0,0,0}};
  for (qp = 0; qp < 16; qp++) {
    POINT param_loc (gps16[qp].r, gps16[qp].s, gps16[qp].t);
    EVALPT evalpt (u, gcell(), param_loc);
    evalpt.eval (geometry);
    int nbfuns = evalpt.nbfuns();
    if ((qp>=0) && (qp<4)) fact = rho * GPW16_1 * evalpt.detJ ();
    else       fact = rho * GPW16_2 * evalpt.detJ ();
    // now loop over all nodes
    for (int I = 0; I < nbfuns; I++) {
      double N_I = evalpt.N (I);
      for (int J = 0; J < nbfuns; J++) {
        double N_J = evalpt.N (J);
        double value =  mmult*fact*N_I*N_J;
        m[0][0] = m[1][1] = m[2][2] = value; //RAW
        evs->add_to_op (dv[J], dv[I], m);
      } 
    }
  }
  return true;
}


static ECELL_MODL *
make (GCELL *gcell)
{
  return (new ECELL_MODL_T10 (gcell));
}

bool
ECELL_MODL_T10::register_make_func ()
{
  return ECELL_MODL::register_make_func (make,
                                                     string (GCELL_SOLID_T10::TYPE_NAME),
                                                     string ("default"));
}

double ECELL_MODL_T10::mass (FIELD_VECTOR *geometry) {
  double mass = 0;
    //CHECK (mass > 0, EXCEPTION_BAD_VALUE,;);
  return mass;
}
