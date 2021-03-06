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
#ifndef PROTO_HEXVUE_H
# define PROTO_HEXVUE_H

#include <list>
#include "algo.h"
#include "gmesh.h"
#include "ecell_hexvue.h"
extern "C" {
#include "ioP.h"
}

template <int NUM_COMPONENTS>
class PROTO_HEXVUE {

 public: // object functions ////////////////////////////////////////

  PROTO_HEXVUE (ALGO *algo);
  ~PROTO_HEXVUE () {
    for (typename list<ECELL_HEXVUE<NUM_COMPONENTS> * >::iterator i = _ecells.begin ();
         i != _ecells.end(); i++) {
      delete *i;
    }
  }

  void write_hexvue_data (ckit_iobuf_t IOB, 
                      FIELD_VECTOR *display_on_geometry,
                      FIELD_VECTOR *eval_on_geometry,
                      FIELD<NUM_COMPONENTS> *field);

 private: // object functions /////////////////////////////////////////

  ECELL_HEXVUE<NUM_COMPONENTS> *make_ecell (GCELL *gcell, string implementation);

 private: // object data /////////////////////////////////////////////

  ALGO                                  *_algo;
  list <ECELL_HEXVUE<NUM_COMPONENTS> * > _ecells;

 private: // object functions ////////////////////////////////////////
  
  void create_ecells (FIELD<NUM_COMPONENTS> *field);

};

template <int NUM_COMPONENTS>
void PROTO_HEXVUE<NUM_COMPONENTS>::write_hexvue_data (ckit_iobuf_t IOB,
                                                      FIELD_VECTOR *display_on_geometry,
                                                      FIELD_VECTOR *eval_on_geometry,
                                                      FIELD<NUM_COMPONENTS> *field)
{
  create_ecells (field);
  typename list<ECELL_HEXVUE<NUM_COMPONENTS> * >::iterator i = _ecells.begin ();
  while (i != _ecells.end ()) {
    (*i)->write (IOB, display_on_geometry, eval_on_geometry, field);
    i++;
  }
}

template <int NUM_COMPONENTS>
void
PROTO_HEXVUE<NUM_COMPONENTS>::create_ecells (FIELD<NUM_COMPONENTS> *field)
{
  SMART_HANDLE<BFUN_SET> bfun_set = field->bfun_set();
  _ecells.clear ();
  // now create the ecells
  BFUN_SET::gsubmesh_enumerator_t sme = bfun_set->gsubmesh_enumerator ();
  GSUBMESH *gsubmesh;
  sme.reset ();
  while ((gsubmesh = sme.next ())) {
    GSUBMESH::gcell_group_enumerator_t gge = gsubmesh->gcell_group_enumerator ();
    GCELL_GROUP *gcell_group;
    gge.reset ();
    while ((gcell_group = gge.next ())) {
      string path = "algorithms/hexvue/"
        + _algo->name () + "/gcell_groups/" + gcell_group->name ();
      string implementation = _algo->mgr()->db()->DB_GET_STRING (path + "/implementation");
      GCELL_GROUP::gcell_enumerator_t ge = gcell_group->gcell_enumerator ();
      GCELL *gcell;
      ge.reset ();
      while ((gcell = ge.next ())) {
        if (gcell->is_topmost_gcell_in_field (field)) { 
          ECELL_HEXVUE<NUM_COMPONENTS> *ecell
            = ECELL_HEXVUE<NUM_COMPONENTS>::make_ecell (gcell, implementation);
          if (ecell)
            _ecells.push_back (ecell);
        }
      }
    }
  }
}


template <int NUM_COMPONENTS>
PROTO_HEXVUE<NUM_COMPONENTS>::PROTO_HEXVUE (ALGO *algo)
{
  // instantiated for (with)
  _algo = algo; CHECK (_algo != NULL, EXCEPTION_NULL_PTR,;);
  _ecells.clear ();
}



#endif
