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
#include "ecell_steady_diffusion.h"
#include "evalpt.h"
#include "select_les.h"
#include "proto_steady_diffusion.h"

OFACTORY<ECELL_STEADY_DIFFUSION, ECELL_STEADY_DIFFUSION::make_ecell_func, GCELL *> * ECELL_STEADY_DIFFUSION::_ecell_factory = 0;

void
ECELL_STEADY_DIFFUSION::attach (class PROTO_STEADY_DIFFUSION *proto, string ggpath)
{
  // Get data from the protocol
  DB *db = proto->db();
  FIELD_SCALAR *phi = proto->phi();
  string phi_name = phi->name();
  EVALPT evalpt = EVALPT (phi, gcell());
  // Bind degrees of freedom to the LES
  evalpt.eval ();
  int nbfuns = evalpt.nbfuns ();
  _dofmappers.resize (nbfuns); 
  LES *les = proto->les();
  les->conn_graph_new_batch();
  for (int i = 0; i < nbfuns; i++) {
    BFUN *bfun = evalpt.bfun (i);
    FEN *fen = bfun->fen ();
    _dofmappers[i] = les->dofmapper (phi_name, phi, phi->dofparam_id (fen));
  }
  les->update_conn_graph();
  // Material
  string matname = db->DB_GET_STRING (ggpath + "/material");
  string mattype = db->DB_GET_STRING ("materials/" + matname + "/type");
  _mat = proto->mat (mattype, matname);
}

bool
ECELL_STEADY_DIFFUSION::register_make_func (ECELL_STEADY_DIFFUSION::make_ecell_func make, string gcell_type, string implementation)
{
  if (_ecell_factory == 0) {
    ECELL_STEADY_DIFFUSION::_ecell_factory = new OFACTORY<ECELL_STEADY_DIFFUSION, ECELL_STEADY_DIFFUSION::make_ecell_func, GCELL* >;
  }
  return (ECELL_STEADY_DIFFUSION::_ecell_factory)->register_make_func (make, gcell_type, implementation);
}


ECELL_STEADY_DIFFUSION *
ECELL_STEADY_DIFFUSION::make_ecell (GCELL *gcell, string implementation)
{
  CHECK (ECELL_STEADY_DIFFUSION::_ecell_factory, EXCEPTION_NULL_PTR,;);
  return (ECELL_STEADY_DIFFUSION::_ecell_factory)->make (gcell, gcell->type_name (), implementation);
}
