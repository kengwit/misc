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
#include "bfun_fe.h"
#include "bfun_set.h"
#include "logger_stream.h" //RAW

BFUN_FE::BFUN_FE (FEN *fen, vector <GCELL *> gcells, BFUN_SET *bfun_set)
  : BFUN (fen,bfun_set,gcells) {}
 


BFUN *BFUN_FE::clone (BFUN_SET *bfun_set)
{
  BFUN_FE *new_bfun = new BFUN_FE (this->fen (), _gcells, bfun_set);  
  return new_bfun;
}


bool 
BFUN_FE::refine (REF_CTX *ref_ctx)
{
  bool result = false;
  for (vector <class GCELL*>::iterator i = _gcells.begin (); i != _gcells.end(); i++) {
    GCELL *gc = *i;
    if (!gc->divided()) { 
      gc->divide(ref_ctx);
      result = true;
    }
  } 
  
  set <FEN*> rf;
  refinement_set(&rf);
//   cout << "*****refinement set of***"  <<this->fen()->id()<<"\t ref loc \t "<<this->fen()->ref_loc()<<"\n";  //RAW:ashi
   for (set <FEN *>::iterator ii = rf.begin(); ii != rf.end(); ii++) {
    BFUN*bfun = bfun_set()->fen_to_bfun_any(*ii);
  //  cout << (*ii)->id() <<"\n"; //RAW:ashi        
    if(bfun==0){ 
      vector<GCELL*> gcells = collect_child_gcells (*ii);
      bfun = new BFUN_FE(*ii,gcells, bfun_set());
    }
    bfun->set_is_active(true);
    bfun_set()->add(bfun);
  }
  if (bfun_set()->refinement_strategy()!=bfun_set()->TRUE_HIERARCHICAL) {
    set_is_active(false);
    bfun_set()->add(this);
  }
 //  cout << "*****end*****"  <<"\n";  //RAW:ashi
  return result;
}


void 
BFUN_FE::unrefine ()
{
  set <FEN*> child_set;
  refinement_set(&child_set);
  if (child_set.size() != 0) { 
    for (set <FEN *>::iterator i = child_set.begin(); i != child_set.end(); i++) {
      BFUN *child_bfun =  this->bfun_set()->fen_to_bfun_any(*i); //RAW 
      if (child_bfun !=0 ) {
        if (!child_bfun->has_inactive_parent(this)) {
          child_bfun->set_is_active(false);
          child_bfun->bfun_set()->add(child_bfun);
        }
      }
    }
  } 
  if (bfun_set()->refinement_strategy()!=bfun_set()->TRUE_HIERARCHICAL) {
    this->set_is_active(true);
    bfun_set()->add(this);
  }
}


bool BFUN_FE::has_inactive_parent(BFUN *current_parent)
{
  
  set<FEN*> parent_set; 
  this->unrefinement_set(&parent_set);
  for (set<FEN*>::iterator i = parent_set.begin();i != parent_set.end(); i++) {
    BFUN *bfun = bfun_set()->fen_to_bfun_any(*i);
    CHECK (bfun!=0, EXCEPTION_NULL_PTR,;);
    if ((!bfun->is_active()) && (bfun != current_parent) ) return true; //RAW
  }
  return false;
}

vector<GCELL*>  BFUN_FE::collect_child_gcells(FEN *fen)
{
  set <GCELL*> result_temp;
  for (sizet k = 0;k<this->ngcells(); k++) {
    GCELL *gcell = this->gcell(k);
    for (sizet i = 0; i<gcell->nchildren(); i++ ){
      GCELL *child_gcell = gcell->child(i);
      for (sizet l=0; l<child_gcell->conn()->nfens(); l++) {
        if (fen == child_gcell->conn()->fen(l)) result_temp.insert(child_gcell);  
      } 
    }
  } 
  vector  <GCELL*> result;
  result.clear();
  for (set<GCELL*>::iterator i = result_temp.begin(); i != result_temp.end();i++) {
    result.push_back(*i);
  }
  return result;  
}
