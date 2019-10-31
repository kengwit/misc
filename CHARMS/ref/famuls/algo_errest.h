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
#ifndef ALGO_ERREST_H
#define ALGO_ERREST_H

#include "fen.h"
#include "algo.h"
#include "gmesh.h"
#include "db.h"
#include "field.h"
#include "mgr.h"
#include "algo_err_rep.h"

/**
   This class defines a trivial interface for error estimation.
*/
template <int NUM_COMPONENTS>
class ALGO_ERREST : public ALGO_ERR_REP {
  
 public: // object methods //@@//@@//@@//@@//@@//@@//@@//@@//@@//@@//
  
  ALGO_ERREST (string name, MGR *mgr) : ALGO_ERR_REP (name, mgr) { }
  virtual ~ALGO_ERREST () {}
  /**
     Assemble `error' for each basis function pair in the field.
     Use geometry given as argument.
   */
  virtual void assemble_error (FIELD_VECTOR *geometry, FIELD<NUM_COMPONENTS> *field) {}
  virtual void assemble_error (FIELD_VECTOR *geometry, FIELD<NUM_COMPONENTS> *field, double ttime) {};
  virtual void accumulate_error (FIELD_VECTOR *geometry, FIELD<NUM_COMPONENTS> *field, double ttime) {};

  /**
    This method returns a scalar field that describes
    the desired mesh size distribution derived from errors collected by assemble error. 
  */
  SMART_HANDLE<FIELD<1> > hbar (SMART_HANDLE<BFUN_SET> bfun_set,
                                sizet target_nbfuns) {
    FIELD<1> *hbar = new FIELD<1>("hbar",bfun_set,0); 
    BFUN_SET::bfun_enumerator_t e = bfun_set->bfun_enumerator ();
    double total_error_vol = 0;
    sizet nbfuns = 0;
    BFUN *bfun;
    e.reset ();
    while ( (bfun = e.next()) != 0) {
      double error = errval (bfun);
      double vol = bfun->char_vol();
      total_error_vol += vol * fabs(error);
      nbfuns++;
    }
    // if target_nbfuns not specified, set it
    if (target_nbfuns == 0) target_nbfuns = nbfuns;
    // desired balanced error density
    double avg_err_per_bfun = (total_error_vol/target_nbfuns);
/*     cerr << endl; */
/*     cerr << "nbfuns " << nbfuns << " target_nbfuns " << target_nbfuns << endl; */
/*     cerr << "avg_err_per_bfun " << avg_err_per_bfun << endl; */
/*     cerr << "total_error_vol " << total_error_vol << endl; */
    init_hbar (hbar, avg_err_per_bfun);
    return hbar; 
  }
    
 private:
  
  void init_hbar (FIELD<1> *h_bar, double avg_err_per_bfun);
  
};

#include "algo_hexvue.h"

  // This is a slight variation on proto field transfer.
template <int NUM_COMPONENTS>
void ALGO_ERREST<NUM_COMPONENTS>::init_hbar (FIELD<1> *h_bar, double avg_err_per_bfun) {
  FIXED_VECTOR<1> z(0);
  GMESH *gmesh = h_bar->bfun_set()->gmesh();
  G3D_box_t box;
  gmesh->box(&box);
  double xr = G3D_BOX_RANGE(box,x);
  double yr = G3D_BOX_RANGE(box,y);
  double zr = G3D_BOX_RANGE(box,z);
  double gmesh_char_dim = max(xr,max(yr,zr));
  
  // Now start the transfer, from the bottom to the top
  sizet maxlevel = h_bar->max_level (h_bar);
  for (sizet level = 0; level < maxlevel+1; level++) {
    const sizet npairs = h_bar->npairs();
    for (sizet j = 0; j < npairs; j++) {
      FIELD_PAIR<1> *fp = h_bar->jth_field_pair (j);
      BFUN *bfun = fp->bfun();
      if (bfun->level() == level) {
        FEN *fen = bfun->fen();
        GCELL *gcell = bfun->gcell(0); // any gcell would do
        // Evaluate basis function values at level lower than the level of bfun
        EVALPT destevalpt (h_bar, gcell, fen);
        destevalpt.eval ();
        double ei = errval(bfun);
        FIXED_VECTOR<1> srcv = 0;
        if (ei != 0) { //This is the formula in use for mesh size
          // RAW this should be the manifold dimension
          double mandim = 3;
          srcv(0) = pow(fabs(avg_err_per_bfun/ei),1./mandim); 
          if (srcv(0) > gmesh_char_dim) srcv(0) = gmesh_char_dim;
        } else {
          srcv(0) = gmesh_char_dim;
        }
        // Loop over all basis functions in the dest field at destevalpt
        double N_j = 1;  // Adjust for PU 
        sizet nbfuns = destevalpt.nbfuns();
        for (sizet j = 0; j < nbfuns; j++) {
          BFUN *dbfun = destevalpt.bfun(j);
          if (dbfun != bfun) {
            double N = destevalpt.N(j);
            BFUN_DOFPARAM_PAIR_ID destdpid = h_bar->dofparam_id (dbfun->fen());
            CHECK (destdpid != INVALID_BFUN_DOFPARAM_PAIR_ID, EXCEPTION_NULL_PTR,;);
            FIXED_VECTOR<1> dp = h_bar->field_pair(destdpid)->dofparam();
            srcv.add (-N, dp);
          } else {
            N_j = destevalpt.N(j);
          }
        }
        srcv.scale(1/N_j); // Adjust for PU 
        fp->set_dofparam_all (srcv);
      } // if
    } // for
  } // for
}

#endif
