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
#include "ref_fen_src.h"

REF_FEN_SRC::REF_FEN_SRC ()
{
  _crlmm.clear ();
}

static sizet smallest_fen_id (CONN_BASE *conn)
{
  sizet fenid = UINT_MAX;
  for (sizet j = 0; j < conn->nfens (); j++) {
    FEN *fen = conn->fen (j);
    if (fen->id () < fenid) fenid = fen->id ();
  }
  return fenid;
}


FEN *REF_FEN_SRC::get_ref_fen (CONN_BASE *conn, sizet indx, FEN_MAKER *fen_maker)
{
  CHECK (conn->nreffens () > 0, EXCEPTION_BAD_ACCESS,;); // why do you want a refinement node for this?

  // First find the min(fenid) -> conn_ref_list map
  conn_ref_list_map_map_t::iterator crlmmi = _crlmm.find (conn->conn_code ());
  if (crlmmi == _crlmm.end ()) {
    _crlmm.insert (conn_ref_list_map_map_t::value_type (conn->conn_code (), new conn_ref_list_map_t()));
    crlmmi = _crlmm.find (conn->conn_code ());
  }
  CHECK (crlmmi != _crlmm.end(), EXCEPTION_NULL_PTR ,;);
  conn_ref_list_map_t *crlm = crlmmi->second;
  // Now we've got the min(fenid) -> conn_ref_list map; find the conn_ref_list
  sizet fen_id_tag = smallest_fen_id (conn);
  conn_ref_list_map_t::iterator crlmi = crlm->find (fen_id_tag);
  if (crlmi == crlm->end ()) {
    crlm->insert (conn_ref_list_map_t::value_type (fen_id_tag, new conn_ref_list_t()));
    crlmi = crlm->find (fen_id_tag);
  }
  CHECK (crlmi != crlm->end (), EXCEPTION_NULL_PTR ,;);
  list <CONN_REF *> *l = crlmi->second;
  // Now we've got the list of conn_ref's
  CONN_REF *conn_ref = 0;
  for (list <CONN_REF *>::iterator li = l->begin (); li != l->end (); li++) {
    conn_ref = *li;
    FEN *fen = conn_ref->get_ref_fen (conn, indx);
    if (fen) return fen;
  }
  // Getting all the way down here means there was no such conn_ref; make it and save it
  conn_ref = new CONN_REF (conn);
  //  cerr << "Adding conn_ref for ";
  //  for (sizet k = 0; k < conn->nfens (); k++) cerr << " " << conn->fen(k)->id();
  //  cerr << endl;
  l->push_back (conn_ref);
  conn_ref->refine (fen_maker);
   FEN *fen = conn_ref->get_ref_fen (conn, indx);
  CHECK (fen, EXCEPTION_NULL_PTR,;);
  return fen;
}

