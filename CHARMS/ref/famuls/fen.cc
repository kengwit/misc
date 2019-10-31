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
#include "fen.h"
#include "gmesh.h"
#include "bfun_set.h"

FEN::FEN (unsigned int id, POINT ref_loc, GMESH *gmesh) : _uniqobj() {
  _id = id;
  _ref_loc = ref_loc;
  _gmesh = gmesh;
  _gmesh->add_fen (this);
 // cerr << "_ref_loc = " << _ref_loc[0] << " " << _ref_loc[1] << " " << _ref_loc[2] << " " << endl;
}
