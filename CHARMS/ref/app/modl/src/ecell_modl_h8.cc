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
#include "ecell_modl_h8.h"
#include "evalpt.h"
#include "square_fixed_matrix.h"
#include "gcell_solid_h8.h"

const char ECELL_MODL_H8::_type_name[] = "modl_h8";

ECELL_MODL_H8::ECELL_MODL_H8 (GCELL *gcell)
  : ECELL_MODL::ECELL_MODL (gcell)
{
}

#include "btcb.h"

bool ECELL_MODL_H8::assemble_stiffness_matrix (FIELD_VECTOR *u,
                                                     EVS_OOOFS *evs, FIELD_VECTOR *geometry)
{
  GAUSS_3d_integration_rule_t *gr = gauss_rule_3d (2);
  vector <DOFMAPPER_VECTOR *> dv = dofmappers ();
  const int nbfuns = dv.size ();
  double K[3][3];
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
        btcb_3d_iso_symm_C (N_xJ, N_yJ, N_zJ, N_xI, N_yI, N_zI, lambda, mu, fact, K);
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

bool ECELL_MODL_H8::assemble_geom_stiffness_matrix (FIELD_VECTOR *u, EVS_OOOFS *evs, FIELD_VECTOR *geometry,double p)
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
  for (qp = 0; qp < gr->npoints; qp++) { // for each integration point
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


bool ECELL_MODL_H8::assemble_geom_stiffness_matrix (FIELD_VECTOR *u, EVS_OOOFS *evs, FIELD_VECTOR *geometry)
{
  //This stiffness matrix assumes that material is homogenous isotropic
  GAUSS_3d_integration_rule_t *gr = gauss_rule_3d (2);
  vector <DOFMAPPER_VECTOR *> dv = dofmappers ();
  const int nbfuns = dv.size ();
  double sigma[3][3];
  for (int ii=0;ii<3; ii++)  for (int jj=0;jj<3; jj++) sigma[ii][jj] = 0;
  double lambda = 0;
  double mu = 0;
  int    qp = 0;
  POINT param_loc (gr->ip[0].xi, gr->ip[0].eta, gr->ip[0].theta);
  POINT at = geometry->evaluate (geometry, this->gcell(), param_loc);
  lambda = mat()->var("lambda", at);
  mu = mat()->var ("mu", at);
  sigma[0][0] = mat()->var ("sigma11", at);
  sigma[1][1] = mat()->var ("sigma22", at);
  sigma[2][2] = mat()->var ("sigma33", at);
  sigma[1][0] = sigma[0][1] = mat()->var ("sigma12", at);
  sigma[2][0] = sigma[0][2] = mat()->var ("sigma13", at);
  sigma[2][1] = sigma[1][2] = mat()->var ("sigma32", at);
  
  for (qp = 0; qp < gr->npoints; qp++) { // for each integration point
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
        double NxI_sigma_NxJ =  (N_xJ*sigma[0][0]+ N_yJ* sigma[0][1]+N_zJ* sigma [0][2])*N_xI +
                                (N_xJ*sigma[1][0]+ N_yJ* sigma[1][1]+N_zJ* sigma [1][2])*N_yI +
	                        (N_xJ*sigma[2][0]+ N_yJ* sigma[2][1]+N_zJ* sigma[2][2] )*N_zI ;
        double k[3][3] = {{NxI_sigma_NxJ*fact, 0,                       0},
                           {0,                 NxI_sigma_NxJ *fact,     0},
                           {0,                 0,                       NxI_sigma_NxJ*fact}};
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
ECELL_MODL_H8::assem_consistent_mass (FIELD_VECTOR *u,EVS_OOOFS *evs, FIELD_VECTOR *geometry, double mmult)
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


double ECELL_MODL_H8::mass (FIELD_VECTOR *geometry) {
  double mass = 0;
  vector <DOFMAPPER_VECTOR *> dv = dofmappers ();
  const int nbfuns = dv.size ();
  GAUSS_3d_integration_rule_t *gr = gauss_rule_3d (2);
  POINT param_loc (gr->ip[0].xi, gr->ip[0].eta, gr->ip[0].theta);
  POINT at = geometry->evaluate (geometry, gcell(), param_loc);
  double rho = mat()->var("rho", at);
  for (int qp = 0; qp < gr->npoints; qp++) {
    POINT param_loc (gr->ip[qp].xi, gr->ip[qp].eta, gr->ip[qp].theta);
    EVALPT evalpt (geometry, gcell(), param_loc);
    evalpt.eval (geometry);
    double fact = rho * gr->ip[qp].Wgt * evalpt.detJ ();
    // now loop over all nodes
    for (int I = 0; I < nbfuns; I++) {
      double N_I = evalpt.N (I);
        mass += fact*N_I;
    }
  }
  return mass;
}

double
ECELL_MODL_H8::strain_energy (FIELD_VECTOR *u, FIELD_VECTOR *geometry)
{
  double sigma[3][3];
  for (int ii=0;ii<3; ii++)  for (int jj=0;jj<3; jj++) sigma[ii][jj] = 0;
  double total_strain_energy  = 0;
  GAUSS_3d_integration_rule_t *gr = gauss_rule_3d (2);
  double detj = 0;
  for (int qp = 0; qp < gr->npoints; qp++) {
    double stress[6]  = {0,0,0,0,0,0};
    double strain[6] = {0,0,0,0,0,0};
    double C_mat[6][6];
    double strain_energy = 0;
    POINT param_loc (gr->ip[qp].xi, gr->ip[qp].eta, gr->ip[qp].theta);
    EVALPT evalpt (u, gcell (), param_loc);
    evalpt.eval (geometry);
    POINT at = geometry->evaluate (geometry, this->gcell(), param_loc);
    double lambda = mat()->var("lambda", at);
    double mu     = mat()->var("mu", at);
    sigma[0][0] = mat()->var ("sigma11", at);
    sigma[1][1] = mat()->var ("sigma22", at);
    sigma[2][2] = mat()->var ("sigma33", at);
    sigma[1][0] = sigma[0][1] = mat()->var ("sigma12", at);
    sigma[2][0] = sigma[0][2] = mat()->var ("sigma13", at);
    sigma[2][1] = sigma[1][2] = mat()->var ("sigma32", at);
    calc_stress (param_loc,u,geometry, stress, strain,lambda,mu);
    mat()->mat_stiffness (param_loc,  C_mat);
    double e_sigma0 =  (sigma[0][0]*strain[0]+ sigma[1][1]*strain[1]+sigma[2][2]*strain[2]
		       +sigma[0][1]*strain[3]+sigma[0][2]*strain[4]+sigma[1][2]*strain[5]);
    double e_sigma  = (stress[0]*strain[0]+ stress[1]*strain[1]+stress[2]*strain[2]
		       +stress[3]*strain[3]+stress[4]*strain[4]+stress[5]*strain[5]);
    double sigma0_dinv_sigma0 = 
    sigma[0][0]*(((lambda+mu)/(mu*(3*lambda+2*mu)))*sigma[0][0]-(lambda/(2*mu*(3*lambda+2*mu)))*(sigma[1][1] + sigma[2][2])) + sigma[1][1]*(((lambda+mu)/(mu*(3*lambda+2*mu)))*sigma[1][1]-(lambda/(2*mu*(3*lambda+2*mu)))*(sigma[0][0] + sigma[2][2])) +  sigma[2][2]*(((lambda+mu)/(mu*(3*lambda+2*mu)))*sigma[2][2]-(lambda/(2*mu*(3*lambda+2*mu)))*(sigma[0][0] + sigma[1][1])) + sigma[0][1]*((1/mu)*sigma[0][1]) + sigma[0][2]*((1/mu)*sigma[0][2]) +  sigma[2][1]*((1/mu)*sigma[2][1]);
    strain_energy = e_sigma0 + 0.5 * e_sigma + 0.5 * sigma0_dinv_sigma0;
    total_strain_energy +=     strain_energy*evalpt.detJ ()* gr->ip[qp].Wgt;
  }
  return total_strain_energy;
}



bool
ECELL_MODL_H8::calc_stress (POINT &param_loc,FIELD_VECTOR *u, FIELD_VECTOR *geometry ,double stress[6], double strain[6],double lambda, double mu)
{

  for (int I = 0; I < 6; I++) { strain[I] = 0; }

  EVALPT evalpt (u, gcell (), param_loc);
  evalpt.eval (geometry);
  const int nbfuns = evalpt.nbfuns();
  // Evaluate derivatives, retrieve displacements, compute strains
  for (int I = 0; I < nbfuns; I++) {
    double N_x = evalpt.N_x (I);
    double N_y = evalpt.N_y (I);
    double N_z = evalpt.N_z (I);
    FIXED_VECTOR<3> u_I = u->field_pair (evalpt.dofparam_id(I))->dofparam ();
    strain[0] += N_x * u_I(0);
    strain[1] += N_y * u_I(1);
    strain[2] += N_z * u_I(2);
    strain[3] += N_y * u_I(0) + N_x * u_I(1);
    strain[4] += N_z * u_I(1) + N_y * u_I(2);
    strain[5] += N_z * u_I(0) + N_x * u_I(2);
  }
  const double l2mu = (lambda + 2*mu);
  stress[0] = l2mu * strain[0] + lambda * (strain[1] + strain[2]);
  stress[1] = l2mu * strain[1] + lambda * (strain[0] + strain[2]);
  stress[2] = l2mu * strain[2] + lambda * (strain[1] + strain[0]);
  stress[3] = mu * strain[3];
  stress[4] = mu * strain[4];
  stress[5] = mu * strain[5];

  return true;
}


static ECELL_MODL *
make (GCELL *gcell)
{
  return (new ECELL_MODL_H8 (gcell));
}

bool
ECELL_MODL_H8::register_make_func ()
{
  return ECELL_MODL::register_make_func (make, string (GCELL_SOLID_H8::TYPE_NAME), string ("default"));
}














