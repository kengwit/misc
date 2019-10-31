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
#ifndef FIELD_H
# define FIELD_H

#include <algorithm>
#include <list>
#include <stack>
#include <vector>
#include <set>
#include <string>
#include "famuls.h"
#include "fen.h"
#include "bfun_dofparam_pair_id.h"
#include "fixed_vector.h"
#include "bfun_set.h"
#include "field_base.h"
#include "field_pair.h"
#include "dofmapper.h"

template <int NUM_COMPONENTS>
class FIELD : public FIELD_BASE {

 public: // object functions ////////////////////////////////////////

  /**
     Make a new field with given name, based on the basis function set,
     initialized to initial_value.
   */
  FIELD (string name, SMART_HANDLE<BFUN_SET> bfun_set,
         FIXED_VECTOR<NUM_COMPONENTS> initial_value)
    : FIELD_BASE (name, bfun_set)
    {
      _pairs.clear ();
      _pairs.resize (bfun_set->nbfuns());
      
      /* build the list of pairs */
      BFUN_SET::bfun_enumerator_t e = bfun_set->bfun_enumerator ();
      BFUN *bfun;
      e.reset ();
      while ((bfun = e.next ())) {
        FIELD_PAIR<NUM_COMPONENTS> *fp = new FIELD_PAIR <NUM_COMPONENTS> (bfun, initial_value);
       // cerr << "Field: adding " << bfun->fen()->id() << " to " << name << endl;
       _pairs[bfun_set->dofparam_id (bfun)] = fp;
      }
      /* and set up the map from gcells to active (ie. non-zero bfun) field pairs */
      build_gcell_to_bfun_dofparam_id_map ();
      // RAW debug printout
      //      for (vector <FIELD_PAIR<NUM_COMPONENTS> *>::iterator i = _pairs.begin ();
      //           i != _pairs.end (); i++) {
      //        FIELD_PAIR<NUM_COMPONENTS> *fp = *i;
      //        if (fp->bfun() == 0) cerr << "Silly bfun == 0" << endl;
      //        else cerr << fp->bfun()->fen()->id() << endl;
      //      }
    }

  /**
     Copy constructor.  ILLEGAL!  Use the copy method instead.
   */
  FIELD (const FIELD & f) {
    EXCEPTION_ILLEGAL_USE e;
    throw e;
  }
  /**
     Assignment constructor.  ILLEGAL!  Use the copy method instead.
  */
  FIELD &operator= (const FIELD & f) {
    EXCEPTION_ILLEGAL_USE e;
    throw e;
  }

  /**
     Copy method.  Creates a new field, which is the exact replica 
     of the field on input.  The new field is a deep copy (in particular,
     the dofparam's have the same value), but it references
     the same BFUN_SET.  The new field is given the name on input.
  */
  FIELD *clone (string name);
  /**
     Set all the degree-of-freedom parameters to the value
     specified on input.  In particular, sending in zero vector
     zeros out all of the field degrees-of-freedom parameters.
   */
  void setto (FIXED_VECTOR<NUM_COMPONENTS> toval)
    {
      for (typename vector <FIELD_PAIR<NUM_COMPONENTS> * >::iterator i = this->_pairs.begin ();
           i != this->_pairs.end (); i++) {
        FIELD_PAIR<NUM_COMPONENTS> *fp = *i;
        fp->set_dofparam (toval);
      }
    }
  /**
     Destructor.
   */
  ~FIELD ()    {
    for (typename vector <FIELD_PAIR<NUM_COMPONENTS> * >::iterator i = _pairs.begin ();
         i != _pairs.end (); i++) {
      FIELD_PAIR<NUM_COMPONENTS> *fp = *i;
      delete fp;
    }
    _pairs.clear();
    for (typename gcell_to_bfun_dofparam_id_map::iterator i = _gcell_to_bfun_dofparam_id_map.begin();
         i != _gcell_to_bfun_dofparam_id_map.end(); i++) {
      (*i).second->clear();
      delete (*i).second;
    }
    _gcell_to_bfun_dofparam_id_map.clear();
  }
  
  /**
     Provide the caller with an opaque identifier for a degree of freedom
     associated in the field with given node.  If this identifier is passed
     to *any* field based on the same basis functions set, the field will be
     able to return the value of the degree of freedom in constant time (i.e. quickly). 
   */
  BFUN_DOFPARAM_PAIR_ID dofparam_id (FEN *fen) {
     BFUN *bfun = bfun_set()->fen_to_bfun(fen);
     if (bfun == 0) return  INVALID_BFUN_DOFPARAM_PAIR_ID;
     else           return  bfun_set()->dofparam_id (bfun); 
  }

  /**
     How many pairs are there in this field?
  */
  sizet npairs () const { return _pairs.size (); }
  /**
     Return a pointer to the field pair associated
     with the identifier on input.   This operation will be done
     efficiently (in constant time).
  */
  FIELD_PAIR<NUM_COMPONENTS> *field_pair (BFUN_DOFPARAM_PAIR_ID dofparam_id) {
    if (dofparam_id == INVALID_BFUN_DOFPARAM_PAIR_ID) { // RAW debug
      cerr << "Field " << this->name() << endl;
      CHECK (dofparam_id != INVALID_BFUN_DOFPARAM_PAIR_ID, EXCEPTION_BAD_ACCESS,;);
    }
    return _pairs[dofparam_id];
  }
  /**
     Return a pointer to the j-th field pair.
     Required: 0<=j<npairs().
  */
  FIELD_PAIR<NUM_COMPONENTS> *jth_field_pair (sizet j) {
    CHECK ((j < _pairs.size()), EXCEPTION_BAD_ACCESS,;);
    return _pairs[j];
  }
  /**
     Return the bfun associated with the identifier on input.   This operation will be done
     efficiently (in constant time).
  */
  BFUN *bfun (BFUN_DOFPARAM_PAIR_ID dofparam_id) const {
    CHECK (dofparam_id != INVALID_BFUN_DOFPARAM_PAIR_ID, EXCEPTION_BAD_ACCESS,;);
    return _pairs[dofparam_id]->bfun();
  }
  /**
     Produce a new dofmapper object.
  */
  DOFMAPPER<NUM_COMPONENTS> *dofmapper (BFUN_DOFPARAM_PAIR_ID dofparam_id) {
    CHECK (dofparam_id != INVALID_BFUN_DOFPARAM_PAIR_ID, EXCEPTION_BAD_ACCESS,;);
    bool constrained[NUM_COMPONENTS];
    _pairs[dofparam_id]->constrained (constrained);
    return new DOFMAPPER<NUM_COMPONENTS> (_pairs[dofparam_id]->bfun()->fen(),
                                          dofparam_id, constrained);
  }
  /**
     Evaluate field at the location of a finite element node.
     This node need not be used by the basis function set of the
     field, which complicates matters somewhat: the field needs
     the figure out to whom to delegate the evaluation of the basis function set.
  */
  FIXED_VECTOR<NUM_COMPONENTS> evaluate (FIELD<3> *geometry, FEN *fen);
  /**
     Evaluate field at the location of a finite element node.
     The node must be used by the gcell given as argument; otherwise the results
     are unpredictable.
  */
  FIXED_VECTOR<NUM_COMPONENTS> evaluate (FIELD<3> *geometry, FEN *fen, GCELL *gcell);
  /**
     Evaluate field at the parametric location <xi> at a gcell.
     This gcell must be a leaf cell (ie. no children).
  */
  FIXED_VECTOR<NUM_COMPONENTS> evaluate (FIELD<3> *geometry, GCELL *gcell, POINT &param_loc);
  FIXED_VECTOR<NUM_COMPONENTS> evaluate (GCELL *gcell, POINT &param_loc) ;//RAW
  /**
     Evaluate field at the node.  The node must be referenced by gcell.
  */
  FIXED_VECTOR<NUM_COMPONENTS> evaluate (FEN *fen, GCELL *gcell);
  /**
     Get the
  */
  PAIRS_ITERATORS pairs_active_over_gcell (GCELL *gcell);
  /**
     Is the field active over cell?
  */
  bool active_over (GCELL *) const;
  /**
     How many free degrees of freedom are there in this field?
  */
  sizet ndofparams () {
    sizet nfree = 0;
    for (typename vector <FIELD_PAIR<NUM_COMPONENTS> *>::iterator i = _pairs.begin ();
         i != _pairs.end (); i++) {
      FIELD_PAIR<NUM_COMPONENTS> *fp = *i;
      if (fp->is_constrained ()) {
        bool c[NUM_COMPONENTS];
        fp->constrained (c);
        for (sizet j = 0; j < NUM_COMPONENTS; j++) if (!c[j]) nfree++;
      } else
        nfree += NUM_COMPONENTS;
    }
    return nfree;
  }    
  /**
     Debug display of a field.
  */
  void debug_display ();
  
  sizet max_level (FIELD<NUM_COMPONENTS> *src) {
    const sizet npairs = src->npairs();
    sizet maxlevel = 0;
    
    for (sizet j = 0; j < npairs; j++) {
      FIELD_PAIR<NUM_COMPONENTS> *fp = src->jth_field_pair (j);
      sizet level = fp->bfun()->level();
      if (level > maxlevel) maxlevel = level;
    }
    return maxlevel;
  }


/*
void  check_consist  ()
{
      for (sizet level = 0; level < max_level(this)+1; level++) {
        const sizet npairs = this->npairs();
        for (sizet j = 0; j < npairs; j++) {
          FIELD_PAIR<NUM_COMPONENTS> *fp = this->jth_field_pair (j);
          BFUN *bfun = fp->bfun();
          if (bfun->level() == level) {
            FEN *fen = bfun->fen();
            GCELL *gcell = bfun->gcell(0); // any gcell would do
            EVALPT geomevalpt (this, gcell, fen);
            geomevalpt.eval ();
            FIXED_VECTOR<NUM_COMPONENTS> srcv(0);
            sizet nbfuns = geomevalpt.nbfuns();
            for (sizet J = 0; J < nbfuns; J++) {
                double N = geomevalpt.N(J);
                BFUN_DOFPARAM_PAIR_ID geomdpid = this->dofparam_id (bfun->fen());
                CHECK (geomdpid != INVALID_BFUN_DOFPARAM_PAIR_ID, EXCEPTION_NULL_PTR,;);
                //calculate summation (n*beta*x), it should be equal to x for linear elements . remove .25 if PU. .25 is partition of 2. -1 is to facilitate substraction 
              srcv.add ((N*bfun->beta()*-1), this->field_pair(geomdpid)->dofparam());
//            srcv.add ((N*-1), this->field_pair(geomdpid)->dofparam()); //without PU           
 }
          //  cout << "n*beta()*x =" << srcv << "\n";
          //  cout << "reference location of evaluation point =" << fen->ref_loc () << "\n";
              srcv.add (fen->ref_loc());
              cout << "fen->ref_loc() - n*beta()*x =" << srcv << "\n";
          } // if
        } // for
      } // for
}
*/

/*
void  check_consist  ()
{
   for (sizet j = 0; j < this->npairs(); j++) {
     FIELD_PAIR<NUM_COMPONENTS> *fp = this->jth_field_pair (j);
     BFUN *bfun = fp->bfun();
     FEN *fen = bfun->fen();
     GCELL *gcell = bfun->gcell(0); // any gcell would do
     EVALPT geomevalpt (this, gcell, fen);
     geomevalpt.eval ();
     FIXED_VECTOR<NUM_COMPONENTS> srcv(0);
     sizet nbfuns = geomevalpt.nbfuns();
     for (sizet J = 0; J < nbfuns; J++) {
       double N = geomevalpt.N(J);
       BFUN *dbfun = geomevalpt.bfun(J);
       BFUN_DOFPARAM_PAIR_ID geomdpid = this->dofparam_id (dbfun->fen());
       CHECK (geomdpid != INVALID_BFUN_DOFPARAM_PAIR_ID, EXCEPTION_NULL_PTR,;);
       //calculate summation (n*beta*x) .  
       srcv.add ((N*-1), this->field_pair(geomdpid)->dofparam());
       // srcv.add ((N*-1), this->field_pair(geomdpid)->dofparam()); //without PU       
     }
     //  cout << "n*beta()*x =" << srcv << "\n";
     //  cout << "reference location of evaluation point =" << fen->ref_loc () << "\n";
     srcv.add (fen->ref_loc());
     cout << "fen->ref_loc() - n*beta()*x =" << srcv << "\n";
   }  
}    
*/


/*
void  check_quad_consist  ()
{
   for (sizet j = 0; j < this->npairs(); j++) {
     FIELD_PAIR<NUM_COMPONENTS> *fp = this->jth_field_pair (j);
     BFUN *bfun = fp->bfun();
     FEN *fen = bfun->fen();
     GCELL *gcell = bfun->gcell(0); // any gcell would do
     EVALPT geomevalpt (this, gcell, fen);
     geomevalpt.eval ();
     FIXED_VECTOR<NUM_COMPONENTS> srcv(0);
     sizet nbfuns = geomevalpt.nbfuns();
     for (sizet J = 0; J < nbfuns; J++) {
       double N = geomevalpt.N(J);
       BFUN *dbfun = geomevalpt.bfun(J);
       BFUN_DOFPARAM_PAIR_ID geomdpid = this->dofparam_id (dbfun->fen());
       CHECK (geomdpid != INVALID_BFUN_DOFPARAM_PAIR_ID, EXCEPTION_NULL_PTR,;);
       //calculate summation (n*beta*x) .
       FIXED_VECTOR<NUM_COMPONENTS>  dp = this->field_pair(geomdpid)->dofparam();
       FIXED_VECTOR<NUM_COMPONENTS> dp_2(dp(0)*dp(0),dp(1)*dp(1),dp(2)*dp(2));   
       srcv.add ((N*-1*0.5), this->field_pair(geomdpid)->dofparam());
       // srcv.add ((N*-1), this->field_pair(geomdpid)->dofparam()); //without PU       
     }
  //   cout << "n*beta()*x*x =" << srcv << "\n";
  //   cout << "reference location of evaluation point =" << fen->ref_loc () << "\n";
       srcv.add (fen->ref_loc());
        cout << "fen->ref_loc() - n*beta()*x*x =" << srcv << "\n";
   }
}

*/


  
 private: // object data //////////////////////////////////////////

  typedef map <GCELL *, set <BFUN_DOFPARAM_PAIR_ID> * > gcell_to_bfun_dofparam_id_map;

  vector <FIELD_PAIR<NUM_COMPONENTS> *> _pairs;
  gcell_to_bfun_dofparam_id_map         _gcell_to_bfun_dofparam_id_map;

 private: // object helpers ///////////////////////////////////////

  void build_gcell_to_bfun_dofparam_id_map ();
  
  class PAIR_COMPARATOR {
  public:
    bool operator() (FIELD_PAIR<NUM_COMPONENTS> *p0, FIELD_PAIR<NUM_COMPONENTS> *p1) const {
      return (p0->bfun()->fen()->id() < p1->bfun()->fen()->id());
    }
  };

};

// Define the type of a scalar field
typedef FIELD<1>     FIELD_SCALAR;
// Define the type of a vector field
typedef FIELD<3>     FIELD_VECTOR;


template <int NUM_COMPONENTS>
void
FIELD<NUM_COMPONENTS>::build_gcell_to_bfun_dofparam_id_map ()
{
#define _M _gcell_to_bfun_dofparam_id_map
  _M.clear ();
  // Add empty lists to the map for all gcells active in the field
  for (sizet j = 0, n = _pairs.size(); j < n; j++) {
    BFUN *bfun = _pairs[j]->bfun ();
    for (sizet k = 0, n = bfun->ngcells(); k < n; k++) {
      gcell_to_bfun_dofparam_id_map::iterator i = _M.find (bfun->gcell (k));
      if (i == _M.end ()) {
        set <BFUN_DOFPARAM_PAIR_ID> * emptyset = new set <BFUN_DOFPARAM_PAIR_ID>; emptyset->clear();
        _M.insert (gcell_to_bfun_dofparam_id_map::value_type (bfun->gcell (k), emptyset));
      }
    }
  }
  // Now for each pair ...
  for (sizet j = 0, n = _pairs.size(); j < n; j++) {
    FIELD_PAIR<NUM_COMPONENTS> *pair = _pairs[j];
    BFUN *bfun = pair->bfun ();
    for (sizet k = 0, n = bfun->ngcells(); k < n; k++) { // ... loop over its cells
      GCELL *gcell = bfun->gcell (k);
      CHECK (_M.find (gcell) != _M.end(), EXCEPTION_BAD_ACCESS,;);
      // First propagate your influence to your ancestors
      GCELL *parent = gcell->parent();
      while (parent) {
        gcell_to_bfun_dofparam_id_map::iterator i = _M.find (parent);
        if (i != _M.end()) { // Parent active in field
          (*i).second->insert (j);
        }
        parent = parent->parent();
      }
      // Now for yourself and your children
      stack <GCELL *> s; 
      s.push (gcell);
      while (!s.empty ()) {
        GCELL *agcell = s.top (); s.pop(); // working on this gcell
        gcell_to_bfun_dofparam_id_map::iterator i = _M.find (agcell);
        if (i != _M.end()) { // This cell is active in the field
          (*i).second->insert (j);
        }
        for (sizet l = 0; l < agcell->nchildren(); l++) s.push (agcell->child(l)); // work on the children too
      }
    } // for each gcell in pair support
  } // for each pair
#undef _M  
}

template <int NUM_COMPONENTS>
typename FIELD<NUM_COMPONENTS>::PAIRS_ITERATORS
FIELD<NUM_COMPONENTS>::pairs_active_over_gcell (GCELL *gcell)
{
  gcell_to_bfun_dofparam_id_map::const_iterator i = _gcell_to_bfun_dofparam_id_map.find (gcell);
  if (i != _gcell_to_bfun_dofparam_id_map.end ()) {
    return PAIRS_ITERATORS((*(i->second)).begin(), (*(i->second)).end());
  } else {
    return PAIRS_ITERATORS(0,0);
  }
}

template <int NUM_COMPONENTS>
bool
FIELD<NUM_COMPONENTS>::active_over (GCELL *gcell) const
{
  gcell_to_bfun_dofparam_id_map::const_iterator i = _gcell_to_bfun_dofparam_id_map.find (gcell);
  if (i != _gcell_to_bfun_dofparam_id_map.end ()) {
    return ( ! (i->second)->empty() );
  } else
    return false;
}

template <int NUM_COMPONENTS>
void
FIELD<NUM_COMPONENTS>::debug_display ()
{
  cerr << "vvv " << this->name() << "##########################################################" << endl;
  cerr << "Field " << this->name() << endl;
  for (sizet j = 0; j < _pairs.size (); j++) {
    FIELD_PAIR<NUM_COMPONENTS> *fp = _pairs[j];
    cerr << j << ": fen " << fp->bfun()->fen()->id() << " ";
    cerr << fp->bfun()->fen()->ref_loc();
    cerr << " dofparam " << fp->dofparam();
    cerr << " dofparam_id " << this->dofparam_id (fp->bfun()->fen());
    if (fp->is_constrained()) {
      bool c[NUM_COMPONENTS]; fp->constrained (c);
      cerr << " constraint[";
      for (int i = 0; i < NUM_COMPONENTS; i++) cerr << " " << c[i];
      cerr << "]";
    }
    cerr << " bfun supported by " << fp->bfun()->ngcells() << " gcells on level " << fp->bfun()->level();
    cerr << endl;
  }
  cerr << "^^^ " << this->name() << "##########################################################" << endl;
}

template <int NUM_COMPONENTS>
FIELD<NUM_COMPONENTS> *
FIELD<NUM_COMPONENTS>::clone (string name)
    {
      FIXED_VECTOR<NUM_COMPONENTS> zero = 0;
      FIELD<NUM_COMPONENTS> *newf = new FIELD<NUM_COMPONENTS> (name, this->bfun_set(), zero);
      for (typename vector <FIELD_PAIR<NUM_COMPONENTS> * >::iterator i = this->_pairs.begin ();
           i != this->_pairs.end (); i++) {
        FIELD_PAIR<NUM_COMPONENTS> *fp = *i;
        BFUN_DOFPARAM_PAIR_ID dpid = this->dofparam_id (fp->bfun()->fen());
        FIELD_PAIR<NUM_COMPONENTS> *newfp = newf->field_pair (dpid);
        newfp->set_dofparam (fp->dofparam ());
      }
      return newf;
    }

#include "evalpt.h"

template <int NUM_COMPONENTS>
FIXED_VECTOR<NUM_COMPONENTS>
FIELD<NUM_COMPONENTS>::evaluate (FIELD_VECTOR *geometry, FEN *fen)
{
  BFUN *bfun = this->bfun_set()->fen_to_bfun (fen);
  if (bfun == 0) { // Node is not used by any bfun in the field bfun_set
    // The trouble here is that unless the fen is associated with a basis function
    // in the field, we don't really know how to get to a gcell that is used by the
    // field. 
    cerr << "Node " << fen->id() << " is not used in field " << this->name () << endl;
    NOT_IMPLEMENTED ("this needs to be done yet");
    FIXED_VECTOR<NUM_COMPONENTS> val = 0;
    return val;
  } else { // Yes, should be able to find a gcell
    CHECK (bfun->ngcells() >= 1, EXCEPTION_BAD_ACCESS,;);
    GCELL *gcell = bfun->gcell(0);
    POINT param_loc;
    gcell->map_fen (fen, &param_loc); // map the fen ...
    POINT mapped_to_param_loc; // ... and find the topmost gcell in field containing this point
    GCELL * evgcell = gcell->map_to_topmost_gcell_in_field (this, param_loc, &mapped_to_param_loc);
    EVALPT evalpt (this, evgcell, mapped_to_param_loc);
    evalpt.eval ();
    const int nbfuns = evalpt.nbfuns (); 
    FIXED_VECTOR<NUM_COMPONENTS> val = 0;
    for (int I = 0; I < nbfuns; I++) {
      double N_I = evalpt.N (I);
      BFUN_DOFPARAM_PAIR_ID dpid = evalpt.dofparam_id (I);
      FIELD_PAIR<NUM_COMPONENTS> *fp = this->field_pair (dpid);
      FIXED_VECTOR<NUM_COMPONENTS> dp = fp->dofparam ();
      val.add (N_I, dp);
    }
    return val;
  }
}

template <int NUM_COMPONENTS>
FIXED_VECTOR<NUM_COMPONENTS>
FIELD<NUM_COMPONENTS>::evaluate (FIELD_VECTOR *geometry, FEN *fen, GCELL *gcell)
{
  POINT param_loc;
  gcell->map_fen (fen, &param_loc); // map the fen ...
  POINT mapped_to_param_loc; // ... and find the topmost gcell in field containing this point
  GCELL * evgcell = gcell->map_to_topmost_gcell_in_field (this, param_loc, &mapped_to_param_loc);
  EVALPT evalpt (this, evgcell, mapped_to_param_loc);
  evalpt.eval ();
  const int nbfuns = evalpt.nbfuns (); 
  FIXED_VECTOR<NUM_COMPONENTS> val = 0;
  for (int I = 0; I < nbfuns; I++) {
    double N_I = evalpt.N (I);
    BFUN_DOFPARAM_PAIR_ID dpid = evalpt.dofparam_id (I);
    FIELD_PAIR<NUM_COMPONENTS> *fp = this->field_pair (dpid);
    FIXED_VECTOR<NUM_COMPONENTS> dp = fp->dofparam ();
    val.add (N_I, dp);
  }
  return val;
}

template <int NUM_COMPONENTS>
FIXED_VECTOR<NUM_COMPONENTS>
FIELD<NUM_COMPONENTS>::evaluate (FIELD_VECTOR *geometry, GCELL *gcell, POINT &param_loc)
{
  POINT mapped_to_param_loc;
  GCELL * evgcell = gcell->map_to_topmost_gcell_in_field (this, param_loc, &mapped_to_param_loc);
  EVALPT evalpt (this, evgcell, mapped_to_param_loc);
  evalpt.eval ();
  const int nbfuns = evalpt.nbfuns (); 
  FIXED_VECTOR<NUM_COMPONENTS> val = 0;
  for (int I = 0; I < nbfuns; I++) {
    double N_I = evalpt.N (I);
    BFUN_DOFPARAM_PAIR_ID dpid = evalpt.dofparam_id (I);
    FIELD_PAIR<NUM_COMPONENTS> *fp = this->field_pair (dpid);
    FIXED_VECTOR<NUM_COMPONENTS> dp = fp->dofparam ();
    val.add (N_I, dp);
  }
  return val;
}

template <int NUM_COMPONENTS>
FIXED_VECTOR<NUM_COMPONENTS>
FIELD<NUM_COMPONENTS>::evaluate (GCELL *gcell, POINT &param_loc) //RAW
{
  POINT mapped_to_param_loc;
  GCELL * evgcell = gcell->map_to_topmost_gcell_in_field (this, param_loc, &mapped_to_param_loc);
  EVALPT evalpt (this, evgcell, mapped_to_param_loc);
  evalpt.eval ();
  const int nbfuns = evalpt.nbfuns (); 
  FIXED_VECTOR<NUM_COMPONENTS> val = 0;
  for (int I = 0; I < nbfuns; I++) {
    double N_I = evalpt.N (I);
    BFUN_DOFPARAM_PAIR_ID dpid = evalpt.dofparam_id (I);
    FIELD_PAIR<NUM_COMPONENTS> *fp = this->field_pair (dpid);
    FIXED_VECTOR<NUM_COMPONENTS> dp = fp->dofparam ();
    val.add (N_I, dp);
  }
  return val;
}

template <int NUM_COMPONENTS>
FIXED_VECTOR<NUM_COMPONENTS>
FIELD<NUM_COMPONENTS>::evaluate (FEN *fen, GCELL *gcell)
{
  POINT param_loc;
  gcell->map_fen (fen, &param_loc); // map the fen ...
  POINT mapped_to_param_loc; // ... and find the topmost gcell in field containing this point
  GCELL * evgcell = gcell->map_to_topmost_gcell_in_field (this, param_loc, &mapped_to_param_loc);
  EVALPT evalpt (this, evgcell, mapped_to_param_loc);
  evalpt.eval ();
  const int nbfuns = evalpt.nbfuns (); 
  FIXED_VECTOR<NUM_COMPONENTS> val = 0;
  for (int I = 0; I < nbfuns; I++) {
    double N_I = evalpt.N (I);
    BFUN_DOFPARAM_PAIR_ID dpid = evalpt.dofparam_id (I);
    FIELD_PAIR<NUM_COMPONENTS> *fp = this->field_pair (dpid);
    FIXED_VECTOR<NUM_COMPONENTS> dp = fp->dofparam ();
    val.add (N_I, dp);
  }
  return val;
}

#endif
