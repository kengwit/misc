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
#include "proto_steady_diffusion.h"

// How to ensure that an algorithm uses the protocol for
// a field defined over the proper set of submeshes, and of the
// proper meaning (scalar field)?
PROTO_STEADY_DIFFUSION::PROTO_STEADY_DIFFUSION (ALGO *algo,// algorithm using this protocol
                                                LES *les, // linear equation solver to use
                                                FIELD_VECTOR *geometry, // defines the geometry
                                                FIELD_SCALAR *phi // scalar diffusion field
                                                ) 
{
  // instantiated for (with)
  _algo = algo; CHECK (_algo != NULL, EXCEPTION_NULL_PTR,;);
  _les  = les; CHECK (_les != NULL, EXCEPTION_NULL_PTR,;);
  _geometry = geometry; CHECK (_geometry != NULL, EXCEPTION_NULL_PTR,;);
  _phi = phi; CHECK (_phi != NULL, EXCEPTION_NULL_PTR,;);
  _mat_map = new MAT_MAP <MAT_DIFFUSION_ISO> (algo->mgr()->mat_mgr());

  // now create the ecells
  BFUN_SET::gsubmesh_enumerator_t sme = _phi->bfun_set()->gsubmesh_enumerator ();
  GSUBMESH *gsubmesh;
  sme.reset ();
  while ((gsubmesh = sme.next ())) {
    // gsubmesh->debug_display ();
    GSUBMESH::gcell_group_enumerator_t gge = gsubmesh->gcell_group_enumerator ();
    GCELL_GROUP *gcell_group;
    gge.reset ();
    while ((gcell_group = gge.next ())) {
      string ggpath = "algorithms/steady_diffusion/"
        + algo->name () + "/gcell_groups/" + gcell_group->name ();
      string implementation = _algo->mgr()->db()->DB_GET_STRING (ggpath + "/implementation");
      GCELL_GROUP::gcell_enumerator_t ge = gcell_group->gcell_enumerator ();
      GCELL *gcell;
      ge.reset ();
      while ((gcell = ge.next ())) {
        if (gcell->is_topmost_gcell_in_field (phi)) { 
          ECELL_STEADY_DIFFUSION *ecell = ECELL_STEADY_DIFFUSION::make_ecell (gcell, implementation);
          CHECK (ecell, EXCEPTION_NULL_PTR,;);
          _ecells.push_back (ecell);
          ecell->attach (this, ggpath);
        }
      }
    }
  }
}

// Here the protocol uses the LES and geometry with which it had been
// created; if a different LES or a different geometry should be used, there
// are analogous routines to assemble things for given LES or geometry.
bool PROTO_STEADY_DIFFUSION::assemble_conductivity_matrix ()
{
  _les->zero_lhs ();
  list<ECELL_STEADY_DIFFUSION *>::iterator i = _ecells.begin ();
  while (i != _ecells.end ()) {
    if (!(*i)->assemble_conductivity_matrix (_phi, _les, _geometry)) return false;
    i++;
  }
  _les->finish_lhs ();
  return true;
}

bool PROTO_STEADY_DIFFUSION::assemble_source_terms ()
{
  _les->zero_rhs ();
  for (list<ECELL_STEADY_DIFFUSION *>::const_iterator ecell = _ecells.begin ();
       ecell != _ecells.end (); ecell++) {
    if (!(*ecell)->assemble_source_terms (_phi, _les, _geometry)) return false;
  }
  _les->finish_rhs ();
  return true;
}

bool PROTO_STEADY_DIFFUSION::solve ()
{
  if (!_les->solve ()) return false;
  return _les->solution (_phi->name(), _phi);
}

MAT_DIFFUSION_ISO *PROTO_STEADY_DIFFUSION::mat (string mattype, string matname)
{
  return _mat_map->mat (mattype, matname);
}

PROTO_STEADY_DIFFUSION::~PROTO_STEADY_DIFFUSION () {
  for (list<ECELL_STEADY_DIFFUSION *>::const_iterator i = _ecells.begin ();
       i != _ecells.end (); i++) {
    ECELL_STEADY_DIFFUSION *ecell = *i;
    delete ecell;
  }
  if (_mat_map) delete _mat_map;
}

