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
#ifndef PROTO_ERREST_H
# define PROTO_ERREST_H

/*
  Error estimation protocol.
 */

#include <list>
#include "algo.h"
#include "gmesh.h"
#include "ecell_errest.h"
#include "mat_mgr.h"

template <int NUM_COMPONENTS>
class PROTO_ERREST {

 public: // object functions ////////////////////////////////////////

  PROTO_ERREST (ALGO *algo);
  ~PROTO_ERREST () {
    for (typename list<ECELL_ERREST<NUM_COMPONENTS> *>::iterator i = _ecells.begin ();
         i != _ecells.end (); i++) {
      delete *i;
    }  
  }
  /**
     Assemble error in the approximation of the field phi
     per basis function.  These numbers may be later retrieved
     by calling error ();
  */
  bool assemble_error (FIELD_VECTOR *geometry, FIELD<NUM_COMPONENTS> *field);
  /**
     Accumulate error for bfun on input.
  */
  void accum (BFUN *bfun, double err, double area);
  /**
     Get error for bfun.
  */
  double errval (BFUN *bfun);
  ALGO *algo () const { return _algo; }
    
 private: // object functions /////////////////////////////////////////

  ECELL_ERREST<NUM_COMPONENTS> *make_ecell (GCELL *gcell, string implementation);
  void make_ecells (FIELD<NUM_COMPONENTS> *field);
  bool do_assemble_error (FIELD_VECTOR *geometry, FIELD<NUM_COMPONENTS> *field);

 private: // object data /////////////////////////////////////////////

  class err_rec_t {
  public:
    err_rec_t (double ierr, double iarea) { err = ierr; area = iarea; }
    double err, area;
  };
  typedef map <BFUN *, err_rec_t> error_map_t;

  ALGO                                    *_algo;
  list <ECELL_ERREST<NUM_COMPONENTS> * >   _ecells;
  error_map_t                              _error_map;
  
};


template <int NUM_COMPONENTS>
PROTO_ERREST<NUM_COMPONENTS>::PROTO_ERREST (ALGO *algo)
{
  _algo = algo; CHECK (_algo != NULL, EXCEPTION_NULL_PTR,;);
  _error_map.clear ();
  _ecells.clear ();
}

template <int NUM_COMPONENTS>
void
PROTO_ERREST<NUM_COMPONENTS>::make_ecells (FIELD<NUM_COMPONENTS> *field)
{
  // now create the ecells
  BFUN_SET::gsubmesh_enumerator_t sme = field->bfun_set()->gsubmesh_enumerator ();
  GSUBMESH *gsubmesh;
  sme.reset ();
  while ((gsubmesh = sme.next ())) {
    GSUBMESH::gcell_group_enumerator_t gge = gsubmesh->gcell_group_enumerator ();
    GCELL_GROUP *gcell_group;
    gge.reset ();
    while ((gcell_group = gge.next ())) {
      string ggpath = "algorithms/errest/"
        + _algo->name () + "/gcell_groups/" + gcell_group->name ();
      string implementation = _algo->mgr()->db()->DB_GET_STRING (ggpath + "/implementation");
      GCELL_GROUP::gcell_enumerator_t ge = gcell_group->gcell_enumerator ();
      GCELL *gcell;
      ge.reset ();
      while ((gcell = ge.next ())) {
        if (gcell->is_topmost_gcell_in_field (field)) { // Leaf cell
          ECELL_ERREST<NUM_COMPONENTS> *ecell
            = ECELL_ERREST<NUM_COMPONENTS>::make_ecell (gcell, implementation);
          if (ecell) {
            _ecells.push_back (ecell); //RAW
             ecell->attach (this, ggpath);
          }
        }
      }
    }
  }
}

template <int NUM_COMPONENTS>
bool PROTO_ERREST<NUM_COMPONENTS>::assemble_error (FIELD_VECTOR *geometry,
                                                   FIELD<NUM_COMPONENTS> *field)
{
  return do_assemble_error (geometry, field);
}

template <int NUM_COMPONENTS>
bool PROTO_ERREST<NUM_COMPONENTS>::do_assemble_error (FIELD_VECTOR *geometry,
                                                      FIELD<NUM_COMPONENTS> *field)
{
  if (_ecells.empty ()) make_ecells (field);
  _error_map.clear();
  // Now fill the map.
  for (typename list<ECELL_ERREST<NUM_COMPONENTS> *>::iterator i = _ecells.begin ();
       i != _ecells.end (); i++) {
    (*i)->assemble_error (field, geometry, this);
  }
  return true;
}

template <int NUM_COMPONENTS>
void PROTO_ERREST<NUM_COMPONENTS>::accum (BFUN *bfun, double err, double area)
{
  typename error_map_t::iterator i = _error_map.find (bfun);
  if (i == _error_map.end ()) {
    _error_map.insert (typename error_map_t::value_type (bfun, err_rec_t(0.0,0.0)));
    i = _error_map.find (bfun);
  }
  //  if (err > i->second) i->second = err;
  i->second.err += err;
  i->second.area += area;
}

template <int NUM_COMPONENTS>
double PROTO_ERREST<NUM_COMPONENTS>::errval (BFUN *bfun)
{
  typename error_map_t::iterator i = _error_map.find (bfun);
  if (i == _error_map.end ()) {
    return 0;
  } else {
    return i->second.err;
    //    if (i->second.area > 0) return i->second.err / i->second.area;
    //    else                    return 0.0;
  }
}

#endif
