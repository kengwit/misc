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
#ifndef ALGO_REFINE_H
#define ALGO_REFINE_H

#include "db.h"
#include "algo.h"
#include "fen.h"
#include "gmesh.h"
#include "ref_ctx.h"
#include "field.h"
#include "algo_errest.h"
#include "logger_stream.h"



class ALGO_REFINE : public ALGO {
  
 public: // object methods //@@//@@//@@//@@//@@//@@//@@//@@//@@//@@//
  
  ALGO_REFINE(string name, MGR *mgr, GMESH *gmesh);
  ~ALGO_REFINE();

   
  /**
    Updates refinement context by  refining basis functions from bottom level of hierarchy 
    to top level of hierarchy and 
    unrefining from top level of hierarchy to bottom level of hierarchy. 

    Description of parameters :
    Adapt method uses following .fpar parameters,

      1. h_over_hbar_ref    : fraction of actual over desired mesh size above which refinement is desired 
      2. h_over_hbar_unref  : fraction of actual over desired mesh size below which unrefinement is desired
      3. ref_fraction       : frcation of no. of basis functions that are considered for refinement.  
               
  */

  template <int NUM_COMPONENTS> 
  SMART_HANDLE <BFUN_SET> adapt(ALGO_ERREST<NUM_COMPONENTS> *ee,
                                  FIELD<NUM_COMPONENTS> *field);

  template <int NUM_COMPONENTS> 
  SMART_HANDLE<BFUN_SET> adapt_quasih(ALGO_ERREST<NUM_COMPONENTS> *ee,
                                        FIELD<NUM_COMPONENTS> *field);

  template <int NUM_COMPONENTS> 
  SMART_HANDLE<BFUN_SET> adapt_trueh(ALGO_ERREST<NUM_COMPONENTS> *ee,
                                        FIELD<NUM_COMPONENTS> *field);
  
  /**
     Build a new basis function set from the basis function set
     given as argument.  This step uses the refinement context
     that needs to have been built prior to this invocation by adapt().
  */
 SMART_HANDLE<BFUN_SET> refine (SMART_HANDLE<BFUN_SET> bfun_set);

 BFUN_SET::REFINEMENT_STRATEGY refinement_strategy () const { 
   if (_true_hierarchical)  return BFUN_SET::TRUE_HIERARCHICAL;
   else                     return BFUN_SET::QUASI_HIERARCHICAL;
 }
  
 private: // object data //@@//@@//@@//@@//@@//@@//@@//@@//@@//@@//@@

  GMESH                  *_gmesh;
  REF_CTX                *_ref_ctx;
  bool                    _true_hierarchical;
  double                  _ref_fraction;
  double                  _h_over_hbar_ref;
  double                  _h_over_hbar_unref;
  sizet                   _max_ref_level;
  sizet                   _target_nbfuns;

 private: // object methods //@@//@@//@@//@@//@@//@@//@@//@@//@@//@@

  class BFUN_ERROR_PAIR {
  public:
    BFUN_ERROR_PAIR () : bfun(0), error(0), didit(false) {}
    BFUN *bfun;
    double error;
    bool didit;
    bool operator<(const BFUN_ERROR_PAIR &a) const {
      return (this->error > a.error);
    }
   }; 
 
  /**
    Check wether a basis function follows refinement/unrefinement rules
  */  
  bool check_rule (int rule,BFUN *bfun);

  /**
    Rule for refinement 
  */  
  bool rule_2 (BFUN *bfun);

  /**
    Rule for unrefinement 
  */
  bool rule_3 (BFUN *bfun);  

};
template <int NUM_COMPONENTS> 
SMART_HANDLE<BFUN_SET> ALGO_REFINE::adapt(ALGO_ERREST<NUM_COMPONENTS> *ee,
                                  FIELD<NUM_COMPONENTS> *field)
{
  if (_true_hierarchical){
    return  adapt_trueh(ee,field); 
  } else {
    return  adapt_quasih(ee,field); 
  } 
}
template <int NUM_COMPONENTS> 
SMART_HANDLE<BFUN_SET> ALGO_REFINE::adapt_trueh(ALGO_ERREST<NUM_COMPONENTS> *ee,
                                  FIELD<NUM_COMPONENTS> *field)
{
  SMART_HANDLE<FIELD<1> > hbar = ee->hbar(field->bfun_set(), _target_nbfuns); 
  SMART_HANDLE<LOGGER_STREAM> lsb = this->mgr()->logger("refine", true);
  SMART_HANDLE<BFUN_SET> new_bfun_set = field->bfun_set()->clone();
  sizet maxlevel = field->max_level(field); 
  vector<list<BFUN*> > buckets(maxlevel+1); 
  BFUN_SET::bfun_enumerator_t e = new_bfun_set->bfun_enumerator ();
  BFUN *bfun;
  sizet bfun_counter = 0; sizet nrefined = 0;
  e.reset ();
  while ((bfun = e.next ()) != 0) {
    BFUN *new_bfun = new_bfun_set->fen_to_bfun_any(bfun->fen());
    buckets[bfun->level()].push_back (new_bfun);
    bfun_counter++;
  } 
  double ntorefine = bfun_counter*_ref_fraction;
  sizet ntoref = 0;
  set <BFUN*> pass2;
  for (list <BFUN*>::iterator i = buckets[0].begin ();
       i != buckets[0].end (); i++) pass2.insert(*i);
  //do refinement
  for (sizet nlevel = 0; nlevel<= maxlevel; nlevel++) {
    for (set<BFUN*>::iterator i = pass2.begin(); i != pass2.end(); i++) {
      BFUN *bfun = *i;
      FIXED_VECTOR<1> hbar_v = hbar->evaluate(bfun->fen(),bfun->gcell(0));
      if ((bfun->char_dim() > hbar_v(0) * _h_over_hbar_ref) && 
          (_max_ref_level > nlevel)) {
        ntoref++;
        bfun->refine (_ref_ctx);
        nrefined ++;
        if (ntorefine > 0 && nrefined >= ntorefine) goto done_refinement;
      }
    }
    pass2.clear();
    if (maxlevel >= nlevel +1) {
      for (list <BFUN*>::iterator i = buckets[nlevel+1].begin ();
           i != buckets[nlevel+1].end (); i++)  pass2.insert(*i);
      for (list <BFUN*>::iterator i = buckets[nlevel].begin ();
           i != buckets[nlevel].end (); i++) {
        BFUN *bfun = *i;
        if (bfun->is_active()) {
          set <FEN*> ref_set;
          bfun->refinement_set(&ref_set);
          for (set<FEN*>::iterator j = ref_set.begin(); j != ref_set.end(); j++) {
            BFUN *bfunx = new_bfun_set->fen_to_bfun_any(*j);
            if (bfunx) {
              set <BFUN*>::iterator i1 = pass2.find(bfunx);
              if ( i1 != pass2.end() ) pass2.erase(i1);         
            }
          } 
        } 
      }
    }
  } 
  
 done_refinement:
  
  *lsb <<"# to refine: "<< ntoref << "; # allowed to be refined: " << ntorefine << "; # refined: " << nrefined << endl_notime;
  
  //do unrefinement
  sizet nunref = 0;
  sizet ntounref = 0;
  
  if (_h_over_hbar_unref <= 1){ // if only unrefinement is required
    for (int nlevel1 = maxlevel; nlevel1>= 0; nlevel1--) {
      set<BFUN*>unref_set; set<BFUN*> ignore_set; 
      for (list <BFUN*>::iterator i1 = buckets[nlevel1].begin ();
           i1 != buckets[nlevel1].end (); i1++) {
        BFUN *bfun = *i1;
        FIXED_VECTOR<1> hbar_v = hbar->evaluate(bfun->fen(),bfun->gcell(0));
        if (bfun->char_dim() <= hbar_v(0) * _h_over_hbar_unref) {
          ntounref++;
          set<FEN*> parent_set;
          bfun->unrefinement_set(&parent_set);
          for (set<FEN*>::iterator i=parent_set.begin(); i!=parent_set.end(); i++) {
            FEN *parent_fen = *i;
            BFUN *parent_bfun = new_bfun_set->fen_to_bfun_any(parent_fen); 
            FIXED_VECTOR<1> hbar_v_parent = hbar->evaluate(parent_bfun->fen(),parent_bfun->gcell(0));
            if ( (hbar_v_parent(0) <= parent_bfun->char_dim()) && 
                 (check_rule(3,parent_bfun)) && (!parent_bfun->is_active())) { 
              unref_set.insert(parent_bfun);
            } else {
              ignore_set.insert(parent_bfun);
            } 
          }
        } else {
          set<FEN*> parent_set;
          bfun->unrefinement_set(&parent_set);
          for (set<FEN*>::iterator i=parent_set.begin(); i!=parent_set.end(); i++) {
            BFUN *parent_bfun = new_bfun_set->fen_to_bfun_any(*i);
            if (parent_bfun != 0) ignore_set.insert(parent_bfun);
          } 
        }  
      }
      for (set<BFUN*>::iterator i=unref_set.begin(); i!=unref_set.end(); i++) {
        if (ignore_set.find(*i) == ignore_set.end()) {
          nunref++; 
          (*i)->unrefine();
        } 
      }
    }
  } 
  
  *lsb << "# to unrefine: " << ntounref << "# unrefined: " << nunref << endl_notime;
  buckets.clear();
  new_bfun_set->commit_clone();
  
  *lsb <<
    "# bfuns on input: " << field->bfun_set()->nbfuns() <<
    "; # bfuns on output: " << new_bfun_set->nbfuns() << endl_notime;
  
  return(new_bfun_set);
  
}

template <int NUM_COMPONENTS> 
SMART_HANDLE<BFUN_SET> ALGO_REFINE::adapt_quasih(ALGO_ERREST<NUM_COMPONENTS> *ee,
                                  FIELD<NUM_COMPONENTS> *field)
{
  SMART_HANDLE<FIELD<1> > hbar = ee->hbar(field->bfun_set(), _target_nbfuns); //RAW :ashi
  bool  nchanges;
  SMART_HANDLE<LOGGER_STREAM> lsb = this->mgr()->logger("refine", true);
  SMART_HANDLE<BFUN_SET> new_bfun_set = field->bfun_set()->clone();
  sizet maxlevel = field->max_level(field); 
  do {
    nchanges = false;
    vector<list<BFUN*> > buckets(maxlevel+1); 
    BFUN_SET::bfun_enumerator_t e = new_bfun_set->bfun_enumerator ();
    sizet bfun_counter = 0; sizet nrefined = 0;
    e.reset ();
    BFUN *bfun = 0;
    while ((bfun = e.next ()) != 0) {
      BFUN *abfun = new_bfun_set->fen_to_bfun_any(bfun->fen());
      buckets[bfun->level()].push_back (abfun);
      bfun_counter++;
    } 
    double ntorefine = bfun_counter*_ref_fraction;
    sizet ntoref = 0;

    set <BFUN*> pass2;
    for (list <BFUN*>::iterator i = buckets[0].begin ();
         i != buckets[0].end (); i++) pass2.insert(*i);
    //do refinement
    
    for (sizet nlevel = 0; nlevel<= maxlevel; nlevel++) {
      for (set<BFUN*>::iterator i = pass2.begin(); i != pass2.end(); i++) {
        BFUN *bfun = *i;
        FEN *fen = bfun->fen();
        sizet fenid = fen->id();
        FIXED_VECTOR<1> hbar_v = hbar->evaluate(bfun->fen(),bfun->gcell(0));
      //  cout << "bfun->fen->id=  "<<bfun->fen()->id()<<" hbar_v=   "<< hbar_v(0) << "char_dim =  "<< bfun->char_dim() << "\n";
        if ((bfun->char_dim() > hbar_v(0) * _h_over_hbar_ref) && 
            (_max_ref_level > nlevel)) {
          ntoref++;
          bfun->refine (_ref_ctx); nchanges = 1;
          new_bfun_set->add(bfun);
          nrefined ++;
          if (ntorefine > 0 && nrefined >= ntorefine) goto done_refinement;
        }
      }
      pass2.clear();
      if (maxlevel >= nlevel +1) {
        for (list <BFUN*>::iterator i = buckets[nlevel+1].begin ();
             i != buckets[nlevel+1].end (); i++)  pass2.insert(*i);
        for (list <BFUN*>::iterator i = buckets[nlevel].begin ();
             i != buckets[nlevel].end (); i++) {
          BFUN *bfun = *i;
          if (bfun->is_active()) {
            set <FEN*> ref_set;
            bfun->refinement_set(&ref_set);
            for (set<FEN*>::iterator j = ref_set.begin(); j != ref_set.end(); j++) {
              BFUN *bfunx = new_bfun_set->fen_to_bfun_any(*j);
              if (bfunx) {
                set <BFUN*>::iterator i1 = pass2.find(bfunx);
                if ( i1 != pass2.end() ) pass2.erase(i1);         
              }
            } 
          } 
        }
      }
    } 
    
  done_refinement:
    
  *lsb <<"# to refine: "<< ntoref << "; # allowed to be refined: " << ntorefine << "; # refined: " << nrefined << endl_notime;
    
    //do unrefinement
    sizet nunref = 0;
    sizet ntounref = 0;
    
    if (_h_over_hbar_unref <= 1) { // unrefinement possible?
      buckets.clear(); buckets.resize(maxlevel + 2);
      BFUN_SET::bfun_enumerator_all_t e1 = new_bfun_set->bfun_enumerator_all ();
      e1.reset ();
      BFUN *bfune;
      while ((bfune = e1.next ()) != 0) { 
        buckets[bfune->level()].push_back (bfune);
      }
      for (int nlevel1 = maxlevel; nlevel1 >= 0; nlevel1--) {  
        set<FEN*> set_to_keep; set_to_keep.clear(); 
        for (list <BFUN*>::iterator i1 = buckets[nlevel1].begin ();
             i1 != buckets[nlevel1].end (); i1++) {  
          BFUN *bfun = *i1;  
          FIXED_VECTOR<1> hbar_v_parent = hbar->evaluate(bfun->fen(),bfun->gcell(0));
          set <FEN *> ref_set; ref_set.clear();
          if ((!bfun->is_active()) &&
              (bfun->char_dim() <= hbar_v_parent(0)*_h_over_hbar_unref ) &&
              bfun->is_refined() 
              ) {
            ntounref++;
            // may unrefine only if either there are no children, or all are active
            bfun->refinement_set(&ref_set);
            bool may_unrefine = true; 
            for (set <FEN *>::iterator i2 = ref_set.begin(); i2 != ref_set.end(); i2++) {
              BFUN *child_bfun = new_bfun_set->fen_to_bfun_any(*i2);
              if (child_bfun) if (!child_bfun->is_active()) { 
                may_unrefine = false;
                break;
              }
            }
            if (may_unrefine) {
              bfun->set_is_active(true);  
              new_bfun_set->add(bfun);
              nunref ++;    
            }  
          }
          // make note of all child functions that would be needed
          if (!bfun->is_active()) { 
            if (ref_set.empty()) bfun->refinement_set(&ref_set);
            for (set <FEN *>::iterator i2 = ref_set.begin(); i2 != ref_set.end(); i2++) set_to_keep.insert(*i2); 
          }
        } 
        for (list <BFUN*>::iterator i3 = buckets[nlevel1+1].begin ();
             i3 != buckets[nlevel1+1].end (); i3++) {
          BFUN *cbfun = *i3;
          if ( (cbfun->is_active()) && (set_to_keep.find(cbfun->fen()) == set_to_keep.end()) ) {
            cbfun->set_is_active(false);
            new_bfun_set->add(cbfun); 
          }      
        }
      }
    }  
    *lsb << "# to unrefine: " << ntounref << "; # unrefined: " << nunref << endl_notime;
    
    buckets.clear();
    BFUN_SET::bfun_enumerator_all_t e2 = new_bfun_set->bfun_enumerator_all ();
    e2.reset (); BFUN *bfunx;
    while ((bfunx = e2.next ()) != 0) { 
      if (bfunx->level()>maxlevel) maxlevel = bfunx->level();
    }
    
  } while (nchanges != 0);
  
  new_bfun_set->commit_clone();

  *lsb <<
    "# bfuns on input: " << field->bfun_set()->nbfuns() <<
    "; # bfuns on output: " << new_bfun_set->nbfuns() << endl_notime;
  
  return(new_bfun_set);
}
  
#endif


