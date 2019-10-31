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
#include "ecell_elasticity.h"
#include "evalpt.h"
#include "proto_elasticity.h"

OFACTORY<ECELL_ELASTICITY, ECELL_ELASTICITY::make_ecell_func, GCELL *> * ECELL_ELASTICITY::_ecell_factory = 0;

void
ECELL_ELASTICITY::attach (PROTO_ELASTICITY *proto, string ggpath)
{
  // Get data from the protocol
  DB *db = proto->db();
  FIELD_VECTOR *u = proto->u();
  string u_name = u->name();
  EVALPT evalpt = EVALPT (u, gcell());
  // Bind degrees of freedom to the LES
  evalpt.eval ();
  int nbfuns = evalpt.nbfuns ();
  _dofmappers.resize (nbfuns); 
  LES *les = proto->les();
  les->conn_graph_new_batch();
  for (int i = 0; i < nbfuns; i++) {
    BFUN *bfun = evalpt.bfun (i);
    FEN *fen = bfun->fen ();
    _dofmappers[i] = les->dofmapper (u_name, u, u->dofparam_id (fen));
  }
  les->update_conn_graph();
  // Material
  string matname = db->DB_GET_STRING (ggpath + "/material");
  string mattype = db->DB_GET_STRING ("materials/" + matname + "/type");
  //cerr << "Requesting " << mattype << endl;
  _mat = proto->mat (mattype, matname);
  // Body load
  if (db->param_defined (ggpath + "/body_load")) {
    string loadname = db->DB_GET_STRING (ggpath + "/body_load");
    string loadtype = db->DB_GET_STRING ("loads/" + loadname + "/type");
    //cerr << "Requesting " << loadtype << endl;
    _body_load = dynamic_cast <LOAD_ON_GCELL_BODY *> (proto->load (loadtype, loadname));
    CHECK (_body_load, EXCEPTION_NULL_PTR,;);
  }
  // Traction load
  if (db->param_defined (ggpath + "/traction_load")) {
    string loadname = db->DB_GET_STRING (ggpath + "/traction_load");
    string loadtype = db->DB_GET_STRING ("loads/" + loadname + "/type");
    //cerr << "Requesting " << loadname << "/" << loadtype << endl;
    _traction_load = dynamic_cast <LOAD_ON_GCELL_TRACTION * > (proto->load (loadtype, loadname));
    CHECK (_traction_load, EXCEPTION_NULL_PTR,;);
  }
}

bool
ECELL_ELASTICITY::register_make_func (ECELL_ELASTICITY::make_ecell_func make, string gcell_type, string implementation)
{
  if (_ecell_factory == 0) {
    ECELL_ELASTICITY::_ecell_factory = new OFACTORY<ECELL_ELASTICITY, ECELL_ELASTICITY::make_ecell_func, GCELL* >;
  }
  return (ECELL_ELASTICITY::_ecell_factory)->register_make_func (make, gcell_type, implementation);
}


ECELL_ELASTICITY *
ECELL_ELASTICITY::make_ecell (GCELL *gcell, string implementation)
{
  CHECK (ECELL_ELASTICITY::_ecell_factory, EXCEPTION_NULL_PTR,;);
  return (ECELL_ELASTICITY::_ecell_factory)->make (gcell, gcell->type_name (), implementation);
}
