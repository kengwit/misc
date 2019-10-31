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
#include "conn_ref.h"
#include "fen_maker.h"

void
CONN_REF::refine (FEN_MAKER *fen_maker)
{
  _ref_fens.clear ();
  sizet nreffens = _conn->nreffens ();
  if (nreffens > 0) {
    _ref_fens.resize(nreffens);
    // Generate refinement nodes.
    for (sizet j = 0; j < nreffens; j++) {
      FEN *fen = (*fen_maker) (j);
      _ref_fens[j] = fen;
    }
  }
}

FEN *CONN_REF::get_ref_fen (CONN_BASE *conn, sizet indx)
{
  return _conn->get_ref_fen (conn, _ref_fens, indx);

}
