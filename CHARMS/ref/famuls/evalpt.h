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
#ifndef EVALPT_H
# define EVALPT_H

#include "gcell.h"
#include "fen.h"
#include "bfun_set.h"
#include "field_base.h"
#include "evalbuf.h"

class EVALPT {

 public: // object functions ////////////////////////////////////////

  /**
     Evaluation point constructor.  Use this constructor when evaluating
     connectivity only.
   */
  EVALPT (FIELD_BASE *field, GCELL *gcell) {
    _field     = field;
    _gcell     = gcell;
    _param_loc = 0;
    _evalbuf   = 0;
    _detJ      = 0;
    _connectivity_only = true;
  }

  /**
     Evaluation point constructor.  Use this constructor when
     evaluating the basis functions at a node.  Since no geometry field
     is available, only the basis function values and their
     derivatives with respect to the parametric coordinates will be computed.
  */
  EVALPT (FIELD_BASE *field, GCELL *gcell, FEN *fen) {
    _field = field;
    _gcell    = gcell;
    _param_loc = 0;
    try {
      _gcell->map_fen (fen, &_param_loc);
    } catch (...) {
      CHECK (0, EXCEPTION_BAD_ACCESS,;);
    }
    _evalbuf  = 0;
    _detJ     = 0;
    _connectivity_only = false;
  }

  /**
     Evaluation point constructor.  The parametric location is supplied
     with respect to the gcell.  Since no geometry field
     is available, only the basis function values and their
     derivatives with respect to the parametric coordinates will be computed.
   */
  EVALPT (FIELD_BASE *field, GCELL *gcell, POINT &param_loc) {
    _field = field;
    _gcell    = gcell;
    _param_loc = param_loc;
    _evalbuf  = 0;
    _detJ     = 0;
    _connectivity_only = false;
  }

  /**
     Destructor.
  */
  ~EVALPT () {
    if (_evalbuf) delete _evalbuf;
  }
  
  /**
     Return the parametric location.
  */
  POINT param_loc () { return _param_loc; }

  /**
     Return the associated bfun_set.
  */
  FIELD_BASE *field () { return _field; }

  /**
     Return the gcell associated with this point.
  */
  GCELL *gcell () { return _gcell; }

  /**
     Evaluate the bfun_set associated with this point.
     Only the basis function values and the derivatives of
     the basis functions with respect to the parametric coordinates
     are going to be evaluated.
  */
  bool eval () {
    if (_evalbuf == 0) {
      _evalbuf = new EVALBUF (_field, _field->pairs_active_over_gcell (_gcell));
    }
    return (_gcell->eval_bfun_set (this) && pu());
  }

  /**
     Evaluate the bfun_set associated with this point.
     The basis function values and the derivatives of
     the basis functions with respect to the spatial coordinates,
     whose parameters are given in the geometry field,
     are going to be evaluated.

     The spatial derivatives of the basis functions are computed from the
     derivatives in the parametric coordinates using the inverse
     Jacobian.  Note 1: This makes sense only if the embedding space
     has the same dimension as the manifold.  Note 2: The derivatives
     in the parametric coordinates are overwritten, and are no longer
     available.
   */
  bool eval (FIELD_VECTOR *geometry) {
    if (_evalbuf == 0) {
      _evalbuf = new EVALBUF (_field, _field->pairs_active_over_gcell (_gcell));
    }
    if (_gcell->eval_bfun_set (this))
      return (pu() && compute_spatial_der (geometry));
    else
      return false;
  }
  
  /**
     Get the number of basis functions after the basis function set
     had been evaluated at the point.
  */
  int nbfuns () { return (*_evalbuf).size(); }
  
  /**
     Get the I-th dofparam id.
  */
  BFUN_DOFPARAM_PAIR_ID dofparam_id (int I) const {
    CHECK_THROW (((I >= 0) && (I < (int) (*_evalbuf).size ())), EXCEPTION_BAD_ACCESS,;);
    return (*_evalbuf)(I).dofparam_id;
  }

  /**
     Get the I-th basis function object.
  */
  BFUN *bfun (int I) const {
    CHECK_THROW (((I >= 0) && (I < (int) (*_evalbuf).size ())), EXCEPTION_BAD_ACCESS,;);
    return (*_evalbuf)(I).bfun;
  }

  FEN *fen (int I) const {
    CHECK_THROW (((I >= 0) && (I < (int) (*_evalbuf).size ())), EXCEPTION_BAD_ACCESS,;);
    return (*_evalbuf)(I).fen;
  }

  /**
     Get the value of the I-th basis function.
  */
  double N (int I) const {
    CHECK_THROW (((I >= 0) && (I < (int) (*_evalbuf).size ())), EXCEPTION_BAD_ACCESS,;);
    return (*_evalbuf)(I).N;
  }

  /**
     Get the derivative of the I-th basis function wrt the J-th
     parametric coordinate.  
  */
  double N_der (int I, int J) const {
    CHECK_THROW (((I >= 0) && (I < (int) (*_evalbuf).size ())), EXCEPTION_BAD_ACCESS,;);
    CHECK_THROW (((J >= 0) && (J < _gcell->manifold_dim())), EXCEPTION_BAD_ACCESS,;);
    return (*_evalbuf)(I).N_der[J];
  }

  /**
     Get the derivative of the I-th basis function wrt the first
     spatial coordinate.
  */
  double N_x (int I) const {
    CHECK_THROW (((I >= 0) && (I < (int) (*_evalbuf).size ())), EXCEPTION_BAD_ACCESS,;);
    return (*_evalbuf)(I).N_der[0];
  }    

  /**
     Get the derivative of the I-th basis function wrt the second
     spatial coordinate.  This function returns zero for 1-manifold cells.
  */
  double N_y (int I) const  {
    CHECK_THROW (((I >= 0) && (I < (int) (*_evalbuf).size ())), EXCEPTION_BAD_ACCESS,;);
    return (*_evalbuf)(I).N_der[1];
  }    

  /**
     Get the derivative of the I-th basis function wrt the third
     spatial coordinate.  This function returns zero for 1- and 2-manifold cells.
  */
  double N_z (int I) const { 
    CHECK_THROW (((I >= 0) && (I < (int) (*_evalbuf).size ())), EXCEPTION_BAD_ACCESS,;);
    return (*_evalbuf)(I).N_der[2];
  }    

  void set_N (int I, double N) {
     (*_evalbuf)(I).N = N;
  }

  void set_N_der (int I, int m, double N_der) {
     (*_evalbuf)(I).N_der[m] = N_der;
  }
  
  /**
     Get the Jacobian at the point.
  */
  double detJ () const { return _detJ; }

  /**
     Should the evaluation look at the connectivity only?  In that
     case, basis function values are not available, only the connectivity.
   */
  bool set_connectivity_only (bool flag = true) {
    bool r = _connectivity_only; _connectivity_only = flag; return r;
  }
  
  /**
     Is this point interested in connectivity only?
  */
  bool connectivity_only () { return _connectivity_only; }


  /**
     Debug display
  */
  void debug_display();
  void debug_display(int line);

 private: // object data ////////////////////////////////////////////

  FIELD_BASE                *_field;
  GCELL                     *_gcell;
  bool                       _connectivity_only;
  POINT                      _param_loc;
  EVALBUF                   *_evalbuf;
  double                     _J1[1][1];
  double                     _inv_J1[1][1];
  double                     _J2[2][2];
  double                     _inv_J2[2][2];
  double                     _J3[3][3];
  double                     _inv_J3[3][3];
  double                     _detJ;

 private: // helper methods ////////////////////////////////////////

  bool compute_spatial_der (FIELD_VECTOR *geometry);
  bool pu ();
  void invert_J (gep_dense_matrix_t &J, gep_dense_matrix_t &invJ);

};

#endif
