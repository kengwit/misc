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
#include "gcell.h"
#include "field.h"
#include "evalpt.h"

bool
EVALPT::compute_spatial_der (FIELD_VECTOR *geometry) {
  // Check that we are doing this in the right space:
  // it doesn't make sense for instance for a 1D manifold embedded
  // in the 3D space.
  const sizet MANIFOLD_DIM = _gcell->manifold_dim();
  gep_dense_matrix_t J; gep_initmat (&J, MANIFOLD_DIM);
  gep_dense_matrix_t inv_J; gep_initmat (&inv_J, MANIFOLD_DIM);
  
  // Initialize the geometry parameters
  const sizet nbfuns = (*_evalbuf).size();
  vector <POINT> x(nbfuns);
  for (sizet i = 0; i < nbfuns; i++) {
    BFUN_DOFPARAM_PAIR_ID dpid = geometry->dofparam_id ((*_evalbuf)(i).fen);
    x[i] = geometry->field_pair (dpid)->dofparam ();
    //cerr << "node " << (*_evalbuf)(i).fen->id() << " x[" << i << "]=" << x[i] << endl;
  }
  
  // Now initialize the Jacobian
  for (sizet k = 0; k < MANIFOLD_DIM; k++)
    for (sizet j = 0; j < MANIFOLD_DIM; j++)
      GEP_ELEM(J,k,j) = 0;
  
  // Compute the Jacobian
  for (sizet i = 0; i < nbfuns; i++) {
    for (sizet k = 0; k < MANIFOLD_DIM; k++) {
      double N_der_k = (*_evalbuf)(i).N_der[k];
      for (sizet j = 0; j < MANIFOLD_DIM; j++)
        GEP_ELEM(J,k,j) += N_der_k * x[i](j);
    }
  }
  
  // Compute the inverse of the Jacobian
  invert_J (J, inv_J);
  //  cerr << "_detJ " << _detJ << endl;
  
  double *N_der = new double [MANIFOLD_DIM];
  for (sizet i = 0; i < nbfuns; i++) {
    for (sizet j = 0; j < MANIFOLD_DIM; j++) { N_der[j] = (*_evalbuf)(i).N_der[j]; }
    for (sizet k = 0; k < MANIFOLD_DIM; k++) {
      double d = 0;
      for (sizet j = 0; j < MANIFOLD_DIM; j++) {
        d += GEP_ELEM(inv_J,k,j) * N_der[j];
      }
      (*_evalbuf)(i).N_der[k] = d; 
    }
  }
  delete [] N_der;
  gep_dispose(&J);
  gep_dispose(&inv_J);
  return true;
}

bool
EVALPT::pu () {
  const sizet nbfuns = (*_evalbuf).size();
  if (nbfuns == 0) return true;
  else if (!((*_evalbuf)(0).bfun->bfun_set()->is_pu())) return true;
  else {
    if (_gcell) {
      const sizet MANIFOLD_DIM = _gcell->manifold_dim();
      for (sizet i = 0; i < nbfuns; i++) {
        double beta = (*_evalbuf)(i).bfun->beta();
        (*_evalbuf)(i).N *= beta;
        for (sizet k = 0; k < MANIFOLD_DIM; k++) {
          (*_evalbuf)(i).N_der[k] *= beta;
        }    
      }
    } else {
      for (sizet i = 0; i < nbfuns; i++) {
        double beta = (*_evalbuf)(i).bfun->beta();
        (*_evalbuf)(i).N *= beta;
      }
    }
  }
  return true;
}

void
EVALPT::invert_J (gep_dense_matrix_t &J, gep_dense_matrix_t &invJ) {
  switch (J.n) {
  case 1:
    _detJ = GEP_ELEM(J,0,0);
    if (_detJ < 0) {
      WARN2 ("Negative detJ", _detJ);
      _detJ = 0;
    } else if (_detJ == 0) {
      WARN ("Zero detJ");
      debug_display(__LINE__);
    }    
    
    if (_detJ == 0) {
      GEP_ELEM(invJ,0,0) = 0;
    } else {
      GEP_ELEM(invJ,0,0) = 1. / _detJ;
    }
    break;
  case 2:
    _detJ = GEP_ELEM(J,0,0) * GEP_ELEM(J,1,1) - GEP_ELEM(J,0,1) * GEP_ELEM(J,1,0);
    
    if (_detJ < 0) {
      WARN2 ("Negative detJ", _detJ);
      _detJ = 0;
    } else if (_detJ == 0) {
      WARN ("Zero detJ");
      debug_display(__LINE__);
    }    
    
    if (_detJ == 0) {
      GEP_ELEM(invJ,0,0) = 0;
      GEP_ELEM(invJ,0,1) = 0;
      GEP_ELEM(invJ,1,0) = 0;
      GEP_ELEM(invJ,1,1) = 0;
    } else {
      double detinvJ  = 1 / _detJ;
      GEP_ELEM(invJ,0,0) =  GEP_ELEM(J,1,1) * detinvJ;
      GEP_ELEM(invJ,0,1) = -GEP_ELEM(J,0,1) * detinvJ;
      GEP_ELEM(invJ,1,0) = -GEP_ELEM(J,1,0) * detinvJ;
      GEP_ELEM(invJ,1,1) =  GEP_ELEM(J,0,0) * detinvJ;
    }
    break;
  case 3:
    _detJ =
      (    (GEP_ELEM(J,0,0) * GEP_ELEM(J,1,1) *GEP_ELEM(J,2,2))
           + (GEP_ELEM(J,0,1) * GEP_ELEM(J,1,2) *GEP_ELEM(J,2,0))
           + (GEP_ELEM(J,0,2) * GEP_ELEM(J,1,0) *GEP_ELEM(J,2,1))
           - (GEP_ELEM(J,0,2) * GEP_ELEM(J,1,1) *GEP_ELEM(J,2,0))
           - (GEP_ELEM(J,0,1) * GEP_ELEM(J,1,0) *GEP_ELEM(J,2,2))
           - (GEP_ELEM(J,0,0) * GEP_ELEM(J,1,2) *GEP_ELEM(J,2,1))
           );
    if (_detJ < 0) {
      debug_display(__LINE__);
      WARN2 ("Negative detJ", _detJ);
      _detJ = 0;
    } else if (_detJ == 0) {
      debug_display(__LINE__);
      WARN ("Zero detJ");
    }    
    
    if (_detJ == 0) {
      GEP_ELEM(invJ,0,0) = 0;
      GEP_ELEM(invJ,0,1) = 0;
      GEP_ELEM(invJ,0,2) = 0;
      GEP_ELEM(invJ,1,0) = 0;
      GEP_ELEM(invJ,1,1) = 0;
      GEP_ELEM(invJ,1,2) = 0;
      GEP_ELEM(invJ,2,0) = 0;
      GEP_ELEM(invJ,2,1) = 0;
      GEP_ELEM(invJ,2,2) = 0;
    } else {
      double detinvJ  = 1. / _detJ;
      
      GEP_ELEM(invJ,0,0) =   detinvJ * (GEP_ELEM(J,1,1)*GEP_ELEM(J,2,2) - GEP_ELEM(J,1,2)*GEP_ELEM(J,2,1)); 
      GEP_ELEM(invJ,0,1) = - detinvJ * (GEP_ELEM(J,0,1)*GEP_ELEM(J,2,2) - GEP_ELEM(J,0,2)*GEP_ELEM(J,2,1)); 
      GEP_ELEM(invJ,0,2) =   detinvJ * (GEP_ELEM(J,0,1)*GEP_ELEM(J,1,2) - GEP_ELEM(J,0,2)*GEP_ELEM(J,1,1)); 
      
      GEP_ELEM(invJ,1,0) = - detinvJ * (GEP_ELEM(J,1,0)*GEP_ELEM(J,2,2) - GEP_ELEM(J,1,2)*GEP_ELEM(J,2,0)); 
      GEP_ELEM(invJ,1,1) =   detinvJ * (GEP_ELEM(J,0,0)*GEP_ELEM(J,2,2) - GEP_ELEM(J,0,2)*GEP_ELEM(J,2,0)); 
      GEP_ELEM(invJ,1,2) = - detinvJ * (GEP_ELEM(J,0,0)*GEP_ELEM(J,1,2) - GEP_ELEM(J,0,2)*GEP_ELEM(J,1,0)); 
      
      GEP_ELEM(invJ,2,0) =   detinvJ * (GEP_ELEM(J,1,0)*GEP_ELEM(J,2,1) - GEP_ELEM(J,1,1)*GEP_ELEM(J,2,0)); 
      GEP_ELEM(invJ,2,1) = - detinvJ * (GEP_ELEM(J,0,0)*GEP_ELEM(J,2,1) - GEP_ELEM(J,0,1)*GEP_ELEM(J,2,0)); 
      GEP_ELEM(invJ,2,2) =   detinvJ * (GEP_ELEM(J,0,0)*GEP_ELEM(J,1,1) - GEP_ELEM(J,0,1)*GEP_ELEM(J,1,0)); 
    }
    break;
  default:
    CHECK(0,EXCEPTION_BAD_VALUE,;);
    break;
  }
}

    
void 
EVALPT::debug_display(int line)
{
  cerr << "Line: " << line << endl;
  debug_display();
}

void 
EVALPT::debug_display()
{
  const sizet nbfuns = (*_evalbuf).size();
  cerr << "Evaluation point" << endl;
  cerr << "nbfuns=" << nbfuns << endl;
  for (sizet j = 0; j < nbfuns; j++) {
    cerr << "FEN " << fen(j)->id() << ": N(" << j << ")=" << N(j) << endl;
  }
  if (_gcell) _gcell->debug_display(); cerr << endl;
  GCELL *parent = _gcell->parent();
  while (parent) {
    parent->debug_display(); cerr << endl;
    parent = parent->parent();
  }
  cerr << "param_loc " << _param_loc(0) << " " << _param_loc(1) << " " << _param_loc(2) << endl;
}
