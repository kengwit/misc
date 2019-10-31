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
#include "ecell_modl.h"
#include "evalpt.h"
#include "evs_ooofs.h"
#include "proto_modl.h"

OFACTORY<ECELL_MODL, ECELL_MODL::make_ecell_func, GCELL *> * ECELL_MODL::_ecell_factory = 0;

void
ECELL_MODL::attach (PROTO_MODL *proto, string ggpath)
{
  // Get data from the protocol
  DB *db = proto->db();
  FIELD_VECTOR *u = proto->u();
  string u_name = u->name();
  EVALPT evalpt = EVALPT (u, gcell());
  // Bind degrees of freedom to the EVS
  evalpt.eval ();
  int nbfuns = evalpt.nbfuns ();
  _dofmappers.resize (nbfuns); 
  EVS_OOOFS *evs = proto->evs();
  evs->conn_graph_new_batch();
  for (int i = 0; i < nbfuns; i++) {
    BFUN *bfun = evalpt.bfun (i);
    FEN *fen = bfun->fen ();
    _dofmappers[i] = evs->dofmapper (u_name, u, u->dofparam_id (fen));
  }
  evs->update_conn_graph();
  // Material
  string matname = db->DB_GET_STRING (ggpath + "/material");
  string mattype = db->DB_GET_STRING ("materials/" + matname + "/type");
  //cerr << "Requesting " << mattype << endl;
  _mat = proto->mat (mattype, matname);
}

bool
ECELL_MODL::register_make_func (ECELL_MODL::make_ecell_func make, string gcell_type, string implementation)
{
  if (_ecell_factory == 0) {
    ECELL_MODL::_ecell_factory = new OFACTORY<ECELL_MODL, ECELL_MODL::make_ecell_func, GCELL* >;
  }
  return (ECELL_MODL::_ecell_factory)->register_make_func (make, gcell_type, implementation);
}


ECELL_MODL *
ECELL_MODL::make_ecell (GCELL *gcell, string implementation)
{
  CHECK (ECELL_MODL::_ecell_factory, EXCEPTION_NULL_PTR,;);
  return (ECELL_MODL::_ecell_factory)->make (gcell, gcell->type_name (), implementation);
}
