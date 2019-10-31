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
#ifndef BFUN_H
# define BFUN_H

#include "fen.h"
#include "gcell.h"
#include "enumerator.h"
 
class BFUN_SET;

class GCELL;

class BFUN {

 public: // object functions ////////////////////////////////////////

  BFUN (FEN *fen, class BFUN_SET *bfun_set,vector <GCELL*> gcells) {
     for(vector<GCELL*>::iterator i = gcells.begin(); i != gcells.end(); i++) _gcells.push_back(*i); 
    _is_active = false;
    _bfun_set = bfun_set;
    _fen      = fen;
    _ref_set.clear();
    _beta = 1; // non-PU
  }

  virtual ~BFUN () {}  
  /**
     Get the finite element node associated with this basis function.
   */
  FEN *fen () { return _fen; }
  
  /**
    Get the bfun_set associated with current basis function.
    Note : Correspondance between basis function set and basis function is one to one. 
  */
  class BFUN_SET *bfun_set () { return _bfun_set; }
  
  /**
     Clone a basis function.
   */
  virtual BFUN *clone (class BFUN_SET *bfun_set)  = 0; 
  
  /**
     Return true if basis function is active.
  */
  bool is_active() const {return _is_active;}
  /**
     set the flag that says wether a function is active or not. 
  */ 
  void set_is_active(bool is_active) {_is_active = is_active;}
  /**
    Return true if basis function was refined.
  */
  bool is_refined ();

  /**
     Return number of gcells covered (at least partially) by the bfun.
  */
  sizet ngcells () const { return _gcells.size (); }
  /**
     Return a gcell.
  */
  class GCELL *gcell (sizet j) const { return _gcells[j]; }
 
 /**
     On which level of the refinement hierarchy is the basis function?
  */
  sizet level () const { return _gcells[0]->level();}
  
  /**
   notify the bfun_set  of functions to become active as a result of the refinement of 
   basis function.  Return true if successful.
  */
  virtual  bool  refine (REF_CTX *ref_ctx) { return true; };
  
  /**
   notify the bfun_set  of functions to become inactive as a result of the unrefinement of 
   basis function.  
   */
  virtual void unrefine () {};
  
  /**
   Returns true if one of the parent of the basis function is inactive.
  */  
  virtual bool has_inactive_parent(BFUN *current_parent) = 0;

  /**
   Returns the characteristic dimension of basis function support.  
  */
  double char_dim();
  /**
   Returns the characteristic volume of basis function support.  
  */
  double char_vol();

  BFUN_DOFPARAM_PAIR_ID dofparam_id();
 
  /**
   return the detail_set or refinement set  of basis function depending on the refinement strtegy of 
   basis function set.
  */ 
   bool  refinement_set (set<FEN*> *s);
   void  complete_refinement_set (set <FEN*> *s) ;
   
  /** 
   return unrefinement_set of a basis function. 
  */
   void  unrefinement_set(set <FEN*> *s) ;

/**
    If a gcell which supports basis function is not divided, basis functions associated
    with that gcell are not added to complete_refinement_set or detail_set.  
   */

  void  detail_set (set <FEN*> *s) ;  

  /**
     Return the coefficient that multiplies this basis function.
     It is one (1) for non-partition-of-unity basis; possibly different
     from one otherwise (i.e. for a partition-of-unity basis).
  */
  double beta () const { return _beta; }
  /**
     Set the coefficient that multiplies this basis function.
     It is one (1) for non-partition-of-unity basis; possibly different
     from one otherwise (i.e. for a partition-of-unity basis).
  */
    void set_beta (double beta);
/*
  void set_beta (double beta) {//RAW ashi
    if (beta == 0) {//RAW ashi
         cout << "beta = " << beta <<"\t fen of bfun=" << _fen->id() << "\t is_active"<<_is_active <<"\n";// RAW ashi
    }//RAW ashi
    CHECK(beta != 0, EXCEPTION_BAD_VALUE,;);
    _beta = beta;
  }

*/
 private: // object data ////////////////////////////////////////////
  
  class BFUN_SET *_bfun_set; 
  FEN            *_fen;
  bool            _is_active;  
  set   <FEN*>    _ref_set;
  double          _beta;

 protected: // object data ////////////////////////////////////////////

  vector<GCELL*>  _gcells;   

};

#endif
