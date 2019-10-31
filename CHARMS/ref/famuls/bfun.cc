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
#include "bfun.h"
#include "bfun_set.h"
#include "ref_ctx.h"
extern "C" {
#include "mathP.h"
}

void BFUN::set_beta (double beta) {//RAW ashi
    if (beta == 0) {//RAW ashi
         cout << "beta = " << beta <<"\t fen of bfun=" << _fen->id() << "\t is_active"<<_is_active <<"\n";// RAW ashi
    }//RAW ashi
//     CHECK(beta != 0, EXCEPTION_BAD_VALUE,;); //RAW:ashi
    _beta = beta;
  }


bool BFUN::is_refined ()
{
  for (sizet j = 0; j < ngcells(); j++) { 
    if (!gcell(j)->divided()) {
      return false;
    }
  } 
  set <FEN*> current_set;
  refinement_set(&current_set);
  if (current_set.size() != 0) {
    for (set<FEN*>::iterator i = current_set.begin(); i != current_set.end(); i++) {
      FEN *fen = *i;
      BFUN *bfun =  bfun_set()->fen_to_bfun_any(fen);
      if (bfun == 0) return false; // one of the children is not part of bfun_set.
      else {
        if (!(bfun->is_active() || bfun->is_refined ())) return false;
      }    
    }
    return true;
  }
  return false;  
}


bool  BFUN::refinement_set (set <FEN*> *s) {
    if (this->bfun_set()->refinement_strategy() == this->bfun_set()->TRUE_HIERARCHICAL ) this->detail_set (s) ; 
    else    this->complete_refinement_set (s)  ;
    return 1;
}


void  BFUN::complete_refinement_set (set <FEN*> *s)
{
  (*s).clear();
  for (vector <class GCELL*>::iterator i = _gcells.begin (); i != _gcells.end(); i++) {
    GCELL *gc = *i; 
    if (gc->divided()) gc->complete_refinement_set (this->fen(),s);
  } 
}   

void  BFUN::detail_set (set <FEN*> *s)
{
  (*s).clear(); 
  for (vector <class GCELL*>::iterator i = _gcells.begin (); i != _gcells.end(); i++) {
    GCELL *gc = *i; 
    if (gc->divided())  gc->detail_set (this->fen(), s);
  }
}

void BFUN::unrefinement_set(set <FEN*> *s)
{
  for (sizet i=0; i<this->ngcells(); i++){
    GCELL *gcell = this->gcell(i);
    GCELL *gcell_parent = gcell->parent();
    if (gcell_parent) { // if gcell do not have parent, it is already bottom  
      for (sizet j=0; j<gcell_parent->conn()->nfens(); j++) {
        FEN *fen = gcell_parent->conn()->fen(j);
        BFUN *bfun = bfun_set()->fen_to_bfun_any(fen);
        if (bfun!=0) { 
          set <FEN*> rf; rf.clear(); //   set<FEN*> parent_ref_set = bfun->refinement_set(); 
           bfun->refinement_set(&rf); 
           if (rf.find(this->fen()) != rf.end()) (*s).insert(fen);
        }       
      }
    }
  }
}



double BFUN::char_dim()
{
  double d = 0;
  for (sizet k=0; k< this->ngcells(); k++) {
    double gcd = this->gcell(k)->char_dim();
    d = max(d, gcd);
  }
  return d; 
}

double BFUN::char_vol()
{
  double v = 0;
  for (sizet k=0; k< this->ngcells(); k++) {
    GCELL *gcell = this->gcell(k);
    double gcd = gcell->char_dim();
    sizet mandim = gcell->manifold_dim();
    v += powi(gcd, mandim);
  }
  return v;
}

BFUN_DOFPARAM_PAIR_ID BFUN::dofparam_id() {
    return _bfun_set->dofparam_id(this);
  }

