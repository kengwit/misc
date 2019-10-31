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
#include "proto_elasticity.h"
#include "mat_elasticity.h"

// How to ensure that an algorithm uses the protocol for
// a field defined over the proper set of submeshes, and of the
// proper meaning (scalar field)?
PROTO_ELASTICITY::PROTO_ELASTICITY (ALGO *algo,// algorithm using this protocol
                                                DB *db, // database
                                                LES *les, // linear equation solver to use
                                                FIELD_VECTOR *geometry, // defines the geometry
                                                FIELD_VECTOR *u // displacement field
                                    ) 
{
  // instantiated for (with)
  _algo = algo; CHECK (_algo != NULL, EXCEPTION_NULL_PTR,;);
  _db = db; CHECK (_db != NULL, EXCEPTION_NULL_PTR,;);
  _les  = les; CHECK (_les != NULL, EXCEPTION_NULL_PTR,;);
  _geometry = geometry; CHECK (_geometry != NULL, EXCEPTION_NULL_PTR,;);
  _u = u; CHECK (_u != NULL, EXCEPTION_NULL_PTR,;);
  _mat_map = new MAT_MAP <MAT_ELASTICITY> (algo->mgr()->mat_mgr());

  // now create the ecells
  BFUN_SET::gsubmesh_enumerator_t sme = _u->bfun_set()->gsubmesh_enumerator ();
  GSUBMESH *gsubmesh;
  sme.reset ();
  while ((gsubmesh = sme.next ())) {
    // gsubmesh->debug_display ();
    GSUBMESH::gcell_group_enumerator_t gge = gsubmesh->gcell_group_enumerator ();
    GCELL_GROUP *gcell_group;
    gge.reset ();
    while ((gcell_group = gge.next ())) {
      string ggpath = "algorithms/elasticity/"
        + algo->name () + "/gcell_groups/" + gcell_group->name ();
      string implementation = db->DB_GET_STRING (ggpath + "/implementation");
      GCELL_GROUP::gcell_enumerator_t ge = gcell_group->gcell_enumerator ();
      GCELL *gcell;
      ge.reset ();
      while ((gcell = ge.next ())) {
        if (gcell->is_topmost_gcell_in_field (u)) { 
          ECELL_ELASTICITY *ecell = ECELL_ELASTICITY::make_ecell (gcell, implementation);
          CHECK (ecell, EXCEPTION_NULL_PTR,;);
          _ecells.push_back (ecell);
          ecell->attach (this, ggpath);
          // cerr << "Made ecell for "; gcell->debug_display (); cerr << endl;
        }
      }
    }
  }
}

// Here the protocol uses the LES and geometry with which it had been
// created; if a different LES or a different geometry should be used, there
// are analogous routines to assemble things for given LES or geometry.
bool PROTO_ELASTICITY::assemble_stiffness_matrix ()
{
  _les->zero_lhs ();
  list<ECELL_ELASTICITY *>::iterator i = _ecells.begin ();
  while (i != _ecells.end ()) {
    if (!(*i)->assemble_stiffness_matrix (_u, _les, _geometry)) return false;
    i++;
  }
  _les->finish_lhs ();
  return true;
}

bool PROTO_ELASTICITY::assemble_load_terms ()
{
  _les->zero_rhs ();
  for (list<ECELL_ELASTICITY *>::const_iterator ecell = _ecells.begin ();
       ecell != _ecells.end (); ecell++) {
    if (!(*ecell)->assemble_prescribed_u_load (_u, _les, _geometry)) return false;
    if (!(*ecell)->assemble_body_load (_u, _les, _geometry)) return false;
    if (!(*ecell)->assemble_tractions_load (_u, _les, _geometry)) return false;
  }
  _les->finish_rhs ();
  return true;
}

bool PROTO_ELASTICITY::solve ()
{
  if (!_les->solve ()) return false;
  return _les->solution (_u->name(), _u);
}


MAT_ELASTICITY *PROTO_ELASTICITY::mat (string mattype, string matname)
{
  return _mat_map->mat (mattype, matname);
}

LOAD_ON_GCELL *PROTO_ELASTICITY::load (string loadtype, string loadname)
{
  return dynamic_cast <LOAD_ON_GCELL *>(this->_algo->mgr()->load_mgr()->load (loadtype, loadname));
}

PROTO_ELASTICITY::~PROTO_ELASTICITY () {
  for (list<ECELL_ELASTICITY *>::const_iterator i = _ecells.begin ();
       i != _ecells.end (); i++) {
    ECELL_ELASTICITY *ecell = *i;
    delete ecell;
  }
  if (_mat_map) delete _mat_map;
}  
