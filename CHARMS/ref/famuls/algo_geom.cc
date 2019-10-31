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
#include "bfun_set.h"
#include "field.h"
#include "algo_geom.h"

ALGO_GEOM::ALGO_GEOM(string name, MGR *mgr) : ALGO (name, mgr) {
}
ALGO_GEOM::~ALGO_GEOM(){
}

void
ALGO_GEOM::init (FIELD_VECTOR *geometry, GMESH *gmesh)
{
  FEN *fen;
  GMESH::fens_enumerator_t e = gmesh->fens_enumerator ();
  e.reset ();
  while ((fen = e.next ())) {
    BFUN_DOFPARAM_PAIR_ID dpid = geometry->dofparam_id (fen);
    if (dpid == INVALID_BFUN_DOFPARAM_PAIR_ID) {
      cerr << "Fen " << fen->id() << " is not used by geometry field" << endl;
      CHECK (dpid != INVALID_BFUN_DOFPARAM_PAIR_ID, EXCEPTION_BAD_VALUE,;);
    }
    FIXED_VECTOR<3> ref_loc = fen->ref_loc ();           // This is correct for the unrefined mesh, but
    geometry->field_pair (dpid)->set_dofparam (ref_loc); // incorrect in general for the refined mesh (non-interpolating)
  }
}

