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
#ifndef BFUN_SET_H
# define BFUN_SET_H

#include <map>
#include <vector>
#include <list>
#include "fen.h"
#include "gcell.h"
#include "bfun_fe.h"
#include "gmesh.h"
#include "gsubmesh.h"
#include "smart_handle.h"



class BFUN_SET {

 public: 

  typedef enum  {TRUE_HIERARCHICAL = 1, QUASI_HIERARCHICAL = 0 } REFINEMENT_STRATEGY;

 public: // creation routines ///////////////////////////////////////

  /**
     Replaces the constructor.  The clients are required to use
     only the smart handle, which makes it possible to share
     the basis functions set without requiring somebody to take
     responsibility of deleting the object when it is no longer
     needed for any field.
  */
  static SMART_HANDLE<BFUN_SET> make (REFINEMENT_STRATEGY refinement_strategy,
                                      list <GSUBMESH *> gsubmeshes) { 
    SMART_HANDLE<BFUN_SET> sh = new BFUN_SET (refinement_strategy,gsubmeshes);

#if defined(USE_PU) && USE_PU
    sh->make_pu(); 
#endif
    return sh;
  }
  
  /**
     Replaces the constructor.  The clients are required to use
     only the smart handle, which makes it possible to share
     the basis functions set without requiring somebody to take
     responsibility of deleting the object when it is no longer
     needed for any field.
  */
  static SMART_HANDLE<BFUN_SET> make (REFINEMENT_STRATEGY refinement_strategy,GMESH *gmesh) {
    SMART_HANDLE<BFUN_SET> sh = new BFUN_SET (refinement_strategy,gmesh);

#if defined(USE_PU) && USE_PU
    sh->make_pu();
#endif
    return sh;
  }
  
 public: // declarations ////////////////////////////////////////////
  
  typedef ENUMERATOR <list <GSUBMESH *>, GSUBMESH> gsubmesh_enumerator_t;
  
  /**
    Needed to ensure that iteration proceeds over active functions only.
  */
  class bfun_enumerator_t: public ENUMERATOR <vector <BFUN *>, BFUN>  {
    
  public:
    bfun_enumerator_t (vector<BFUN*> bfuns) : ENUMERATOR <vector <BFUN *>, BFUN>(bfuns) { }
    BFUN *next () {
      BFUN *bfun = 0;
      do {
        bfun = ENUMERATOR <vector <BFUN *>, BFUN>::next ();
      } while ((bfun != 0) && (!bfun->is_active ()));
      return (bfun);
    }
  };


  class bfun_enumerator_all_t: public ENUMERATOR <vector <BFUN *>, BFUN>  {
    
  public:
    bfun_enumerator_all_t (vector<BFUN*> bfuns) : ENUMERATOR <vector <BFUN *>, BFUN>(bfuns) { }
    BFUN *next () {
      return (ENUMERATOR <vector <BFUN *>, BFUN>::next ());
    }
  }; 
  
  
 public: // object functions ////////////////////////////////////////
  
   void debug_display(char filename[20]);
   void debug_display();
  // void check_consist  (FIELD<NUM_COMPONENTS> *geometry);
  
  /**
     Public destructor.
   */
  ~BFUN_SET ();
 
  /**
     Map a finite element node to a active basis function; If there 
     is no active functions in basis function set,0 is returned.

     Remark:  
     This function presents a problem: what if it is desired that
     a bfun_set have several bfuns associated with a single node.
     Does that make sense?  Maybe, but in any case, this dilemma
     may be avoided if each basis function refers to a single node.
  */
  BFUN *fen_to_bfun (FEN *fen);
  /*
  Map finite element node to any basis function in the set, active or inactive; 
  zero is returned if there is no such function.  
  */
  BFUN *fen_to_bfun_any (FEN *fen);
  /**
     Get the number of active basis functions in this set.
  */
  sizet nbfuns () {
    return _nbfuns_active;
  } 
 
  
  /**
     Enumerator of submeshes.
  */
  gsubmesh_enumerator_t gsubmesh_enumerator () {
    gsubmesh_enumerator_t e (_gsubmeshes);
    return e;
  }

  /**
        Enumerator of active basis functions.
   */
  bfun_enumerator_t bfun_enumerator() {
    bfun_enumerator_t e = bfun_enumerator_t (_bfuns);
    return e;
  }
  bfun_enumerator_all_t bfun_enumerator_all() {
    bfun_enumerator_all_t e = bfun_enumerator_all_t (_bfuns);
    return e;
  }
  /**
     Debugging aid: print out the bfun set.
   */
  void display ();
  
  /**
     Get the bfun/dofparam id.  It is the (opaque) id that fields
     use to access pairs.  All fields based on the same bfun_set
     will have the same order of pairs, and hence all may use the
     dofparam identifier generated by the bfun set. INVALID_DOFPARAM_ID
     is returned if basis function is not active.

     Remark : Bfun has to be associated with given basis function set.
  */
  BFUN_DOFPARAM_PAIR_ID dofparam_id (BFUN *bfun);
  /**
     Return the gmesh with which this bfun_set associates.
  */
  GMESH *gmesh() const { return _gmesh; }
 
  /**
   Return one if refinement strategy is true hierarchical. Otherwise it returns 0.
  */
  REFINEMENT_STRATEGY refinement_strategy () const { return _refinement_strategy;}

  /**
    Returns true if basis function passed as argument exists as inactive basis function in 
    basis function set.
  */
  bool in_inactive_set (BFUN *bfun) { 
    return ( (this->dofparam_id(bfun->fen())>= 0) && (!in_active_set(bfun)) );  
  }
  
  /**
    Returns true if basis function passed as argument exists as active basis function in 
    basis function set.
  */
  bool in_active_set (BFUN *bfun) {
    BFUN_DOFPARAM_PAIR_ID dpid = this->dofparam_id(bfun->fen());
    return ( ( dpid >= 0) && (dpid < (BFUN_DOFPARAM_PAIR_ID) _nbfuns_active ));
  }
  
  /**
  make a clone of the basis function set
  */
  SMART_HANDLE <BFUN_SET> clone ();

  /**
  add a basis function to basis function set. This function can be used to change the attribute 
  of a basis function from active to inactive or vice versa.  
  */
  void add (BFUN *bfun);
  
  void commit_clone();

  /**
     Is this a partition-of-unity basis function set?
  */
  bool is_pu() const { 
#if defined(USE_PU) && USE_PU
    return _is_pu; 
#else
    return false;
#endif
  }

#if defined(USE_PU) && USE_PU
  void make_pu();
#endif

 private: // object data ////////////////////////////////////////////

#if defined(USE_PU) && USE_PU
  bool                            _is_pu;
#endif
  GMESH                          *_gmesh;
  list <GSUBMESH *>               _gsubmeshes;
  vector <BFUN *>                 _bfuns;
  vector <BFUN_DOFPARAM_PAIR_ID>  _bfun_dofparam_ids;
  REFINEMENT_STRATEGY             _refinement_strategy;
  bool                            _in_constructor;
  sizet                           _nbfuns_active;
  sizet                           _nbfuns_total;
  sizet                           _mem_raise_size ;
   
 private: // object functions ///////////////////////////////////////

  
  bool is_present(BFUN*bfun);
  void insert_bfun (BFUN *bfun);
  void swap(sizet temp_loc1,sizet temp_loc2);
  void build_from_gcells (list <GCELL *> igcells, BFUN_SET *bfun_set);
  list <GCELL *> build_gcell_list ();
  void build_bfun_dofparam_ids ();
 
  
  

 private: // private constructors ////////////////////////////////////

  /**
     Default: not useful.
  */
  BFUN_SET (REFINEMENT_STRATEGY refinement_strategy);
  /**
     Define basis functions set on the domain covered by the submeshes.
     This constructor takes a list of submeshes and construct FE
     basis functions on all the gcells included in those submeshes.
     The refinement context may be passed in as null (0).  In that case
     all nodes from those submeshes are assumed to be active.
     If the refinement context is passed in as non-null, it is queried
     to find out which nodes should be active.
  */
  BFUN_SET (REFINEMENT_STRATEGY refinement_strategy,list <GSUBMESH *> gsubmeshes);
  /**
     Define a basis function set on the whole geometric mesh.
     The refinement context may be passed in as null (0).  In that case
     all nodes from those submeshes are assumed to be active.
     If the refinement context is passed in as non-null, it is queried
     to find out which nodes should be active.
  */
  BFUN_SET (REFINEMENT_STRATEGY refinement_strategy,GMESH *gmesh);
  BFUN_DOFPARAM_PAIR_ID dofparam_id (FEN *fen);
  
};





#endif
