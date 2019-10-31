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
#include "ecell_steady_diffusion_l2.h"
#include "gcell_line_l2.h"
#include "evalpt.h"
#include "square_fixed_matrix.h"
#include "gaussP.h"

const char ECELL_STEADY_DIFFUSION_L2::TYPE_NAME[] = "steady_diffusion_l2";

ECELL_STEADY_DIFFUSION_L2::ECELL_STEADY_DIFFUSION_L2 (GCELL *gcell) :
  ECELL_STEADY_DIFFUSION::ECELL_STEADY_DIFFUSION (gcell) 
{
}

bool
ECELL_STEADY_DIFFUSION_L2::assemble_conductivity_matrix (FIELD_SCALAR *phi,
                                                         LES *les,
                                                         FIELD_VECTOR *geometry)
{
  GAUSS_1d_integration_rule_t *gr = gauss_rule_1d (1);
  vector <DOFMAPPER_SCALAR *> dv = dofmappers ();
  const int nbfuns = dv.size ();
  SQUARE_FIXED_MATRIX k (nbfuns);
  // First, compute the conductivity matrix
  k.zero ();
  for (int qp = 0; qp < gr->npoints; qp++) {
    POINT param_loc (gr->ip[qp].xi, 0, 0); 
    EVALPT evalpt (phi, gcell (), param_loc);
    evalpt.eval (geometry);
    int nbfuns = evalpt.nbfuns (); 
    double sum_N = 0; //RAW
    //for (sizet i= 0; i<nbfuns; i++) { //RAW
    //  cout <<"N = "<<evalpt.N(i)<<"\t";
    //  sum_N = sum_N + evalpt.N(i); //RAW 
    // } //RAW
    // if (fabs(sum_N)-1 >  1e-7) { //RAW
    //  cout << "*********tolerance = \t"<< fabs(sum_N)-1<<"\t nbfuns\t"<<nbfuns<<"\n"; //RAW
      //CHECK(fabs(sum_N)-1 < 1e-7,EXCEPTION_BAD_VALUE,;); //RAW  
    // } 

    CHECK (evalpt.nbfuns () == nbfuns, EXCEPTION_BAD_VALUE,;);
    POINT at = geometry->evaluate (geometry, this->gcell(), param_loc);
    double fact = mat()->conductivity(at) * gr->ip[qp].Wgt * evalpt.detJ (); 
    // now loop over all nodes 
    for (int I = 0; I < nbfuns; I++) {
      double N_xI = evalpt.N_x (I);
      for (int J = 0; J < nbfuns; J++) {
        double N_xJ = evalpt.N_x (J);
        k(I,J) += (N_xI * N_xJ) * fact;
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
  // and clean up
  gauss_free_1d (gr);
  return true;
}

bool
ECELL_STEADY_DIFFUSION_L2::assemble_source_terms (FIELD_SCALAR *phi,
                                                  LES *les,
                                                  FIELD_VECTOR *geometry)
{
  GAUSS_1d_integration_rule_t *gr = gauss_rule_1d (1);
  vector <DOFMAPPER_SCALAR *> dv = dofmappers ();
  const int nbfuns = dv.size ();
  SQUARE_FIXED_MATRIX k (nbfuns);
  // cerr << "Integrating source terms over " << this->gcell()->conn()->fen(0)->id() << " " << this->gcell()->conn()->fen(1)->id() << endl;
  // First, compute the conductivity matrix
  k.zero ();
  for (int qp = 0; qp < gr->npoints; qp++) {
    POINT param_loc (gr->ip[qp].xi, 0, 0);
    POINT at = geometry->evaluate (geometry, this->gcell(), param_loc);
    EVALPT evalpt (phi, gcell (), param_loc);
    evalpt.eval (geometry);
    int nbfuns = evalpt.nbfuns (); 
    CHECK (evalpt.nbfuns () == nbfuns, EXCEPTION_BAD_VALUE,;);
    // now loop over all nodes 
    for (int I = 0; I < nbfuns; I++) {
      bool constrained[1];
      dv[I]->constrained (constrained);
      if (constrained[0]) {       // need to include the contribution to the rhs
        double N_xI = evalpt.N_x (I);
        BFUN_DOFPARAM_PAIR_ID dpid = dv[I]->dofparam_id ();
        FIXED_VECTOR<1> prescribed_phi = phi->field_pair (dpid)->dofparam ();
        // cerr << "prescribed_phi " << prescribed_phi[0] << endl;
        double fact = (-1) * mat()->conductivity(at) * prescribed_phi[0] * gr->ip[qp].Wgt * evalpt.detJ ();
        for (int J = 0; J < nbfuns; J++) {
          double N_xJ = evalpt.N_x (J);
          double values[1];
          values[0] = (N_xI * N_xJ) * fact;
          les->rhs_add (dv[J], values);
        }
      }
    }
    // Handle internal heat generation
    double fact = mat()->internal_heat_generation_density(at) * gr->ip[qp].Wgt * evalpt.detJ ();
    // cerr << fact << " " << _internal_heat_generation_density << " " << gr->ip[qp].Wgt << " " << evalpt.detJ () << endl;
    for (int I = 0; I < nbfuns; I++) {
      double N = evalpt.N (I);
      double values[1];
      values[0] = N * fact;
      // cerr << "fen " << evalpt.bfun(I)->fen()->id() << " N = " << N << " values[0] = " << values[0] << endl;
      les->rhs_add (dv[I], values);
    }
  }
  // and clean up
  gauss_free_1d (gr);
  return true;
}


static ECELL_STEADY_DIFFUSION *
make (GCELL *gcell)
{
  return (new ECELL_STEADY_DIFFUSION_L2 (gcell));
}

bool
ECELL_STEADY_DIFFUSION_L2::register_make_func ()
{
  return ECELL_STEADY_DIFFUSION::register_make_func (make, string (GCELL_LINE_L2::TYPE_NAME), string ("default"));
}
