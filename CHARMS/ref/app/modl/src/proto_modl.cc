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
#include "proto_modl.h"
#include "mat_elasticity.h"

// How to ensure that an algorithm uses the protocol for
// a field defined over the proper set of submeshes, and of the
// proper meaning (scalar field)?
PROTO_MODL::PROTO_MODL (ALGO *algo,// algorithm using this protocol
                                                DB *db, // database
                                                EVS_OOOFS *evs, // linear equation solver to use
                                                FIELD_VECTOR *geometry, // defines the geometry
                                                FIELD_VECTOR *u // displacement field
                                    ) 
{
  // instantiated for (with)
  _algo = algo; CHECK (_algo != NULL, EXCEPTION_NULL_PTR,;);
  _db = db; CHECK (_db != NULL, EXCEPTION_NULL_PTR,;);
  _evs  = evs; CHECK (_evs != NULL, EXCEPTION_NULL_PTR,;);
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
      string ggpath = "algorithms/modl/"
        + algo->name () + "/gcell_groups/" + gcell_group->name ();
      string implementation = db->DB_GET_STRING (ggpath + "/implementation");
      GCELL_GROUP::gcell_enumerator_t ge = gcell_group->gcell_enumerator ();
      GCELL *gcell;
      ge.reset ();
      while ((gcell = ge.next ())) {
        if (gcell->is_topmost_gcell_in_field (u)) { 
          ECELL_MODL *ecell = ECELL_MODL::make_ecell (gcell, implementation);
          CHECK (ecell, EXCEPTION_NULL_PTR,;);
          _ecells.push_back (ecell);
          ecell->attach (this, ggpath);
          // cerr << "Made ecell for "; gcell->debug_display (); cerr << endl;
        }
      }
    }
  }
  _evs->zero_k();
  _evs->zero_m();
}

// Here the protocol uses the EVS and geometry with which it had been
// created; if a different EVS or a different geometry should be used, there
// are analogous routines to assemble things for given EVS or geometry.
bool PROTO_MODL::assemble_stiffness_matrix (double pressure, double shift)
{
  _evs->assemble_to_K();
  list<ECELL_MODL *>::iterator i = _ecells.begin ();
  while (i != _ecells.end ()) {
    if (!(*i)->assemble_stiffness_matrix (_u, _evs, _geometry)) return false;
    if (!(*i)-> assemble_geom_stiffness_matrix (_u,_evs,_geometry, pressure)) return false;
    if (shift != 0)
      if (!(*i)->assem_consistent_mass(_u, _evs, _geometry, -shift)) return false;
    i++;
  }
   _evs->finish_k ();
   // _evs->dump_lhs ( "stiffness.txt"); 
   return true;
}

bool PROTO_MODL::assemble_stiffness_matrix (double shift)
{
  _evs->assemble_to_K();
  list<ECELL_MODL *>::iterator i = _ecells.begin ();
  while (i != _ecells.end ()) {
    if (!(*i)->assemble_stiffness_matrix (_u, _evs, _geometry)) return false;
    if (!(*i)-> assemble_geom_stiffness_matrix (_u,_evs,_geometry)) return false;
    if (shift != 0)
      if (!(*i)->assem_consistent_mass(_u, _evs, _geometry, -shift)) return false;
    i++;
  }
   _evs->finish_k ();
   // _evs->dump_lhs ( "stiffness.txt"); 
   return true;
}



bool PROTO_MODL::assemble_mass_matrix ()
{
  _evs->assemble_to_M();
  //  SMART_HANDLE<LOGGER_STREAM> lsb = this->_algo_modl->mgr()->logger("assemble consistent mass", true);
  
  //double nz_per_row = 0.01;
   //fles_solver_ctx_t   slesctx = fles_new_spetsc_solver (this, _tot_of_equations, _error_notify);
  
  //slesctx->zero_lhs (slesctx);
  _evs->zero_m ();
  for (list<ECELL_MODL *>::const_iterator ecell = _ecells.begin ();
       ecell != _ecells.end (); ecell++) {
    if (!(*ecell)->assem_consistent_mass(_u, _evs, _geometry, 1.0)) return false;
  }
  _evs->finish_m ();
  // _evs->dump_lhs ( "mass.txt"); 
  //  *lsb << "done"  << endl;
  return true;
}



MAT_ELASTICITY *PROTO_MODL::mat (string mattype, string matname)
{
  return _mat_map->mat (mattype, matname);
}

PROTO_MODL::~PROTO_MODL () {
  for (list<ECELL_MODL *>::const_iterator i = _ecells.begin ();
       i != _ecells.end (); i++) {
    ECELL_MODL *ecell = *i;
    delete ecell;
  }
  if (_mat_map) delete _mat_map;
}  


bool 
PROTO_MODL::solve(double shift)
{
  return(_evs->solve(shift));
}

double PROTO_MODL::total_strain_energy_in_gcell_group (FIELD_VECTOR *u, GCELL_GROUP *gg) {
  double strain_energy = 0;
  for (list <ECELL_MODL *>:: iterator i = _ecells.begin(); i != _ecells.end(); i++) {
    if ((*i)->gcell()->gcell_group() == gg) {  
      strain_energy += (*i)->strain_energy(u, _geometry);
    }
  }
  return strain_energy;
}

double PROTO_MODL::total_mass() {
    double total_mass = 0;
    for (list <ECELL_MODL *>:: iterator i = _ecells.begin(); i != _ecells.end(); i++) {
    	total_mass += (*i)->mass(_geometry);
    }
    return total_mass;
}









