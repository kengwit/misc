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
#include "gcell_point_p1.h"

const char GCELL_POINT_P1::TYPE_NAME[] = "point_p1";

GCELL_POINT_P1::GCELL_POINT_P1 (vector <FEN *> fens) : GCELL ()
{
  sizet nfens = _conn.nfens ();
  CHECK (fens.size () == nfens, EXCEPTION_BAD_VALUE ,;);
  for (unsigned int i = 0; i < nfens; i++) _conn.fen(i) = fens[i];
}

static GCELL *
make (vector <FEN *> fens)
{
  return (new GCELL_POINT_P1 (fens));
}

bool GCELL_POINT_P1::register_make_func ()
{
  return GCELL::register_make_func (make, GCELL_POINT_P1::type_name (), string ("default"));
}

bool
GCELL_POINT_P1::eval_bfun_set (EVALPT *evalpt)
{
  return true;                  // nothing to be done
}

void
GCELL_POINT_P1::map_fen (FEN *fen, POINT *param_loc)
{
  if (fen->id() == _conn.fen(0)->id()) {
    *param_loc = 0;
    return;
  }
  throw EXCEPTION_BAD_VALUE();  // the node is unknown to this cell; raise hell
}

bool
GCELL_POINT_P1::map_to_child (POINT &param_loc, GCELL **child, POINT *child_param_loc)
{
  NOT_IMPLEMENTED (hi);
  return false; // no children
}


bool
GCELL_POINT_P1::map_to_parent (POINT &param_loc, GCELL **parent, POINT *parent_param_loc)
{
  NOT_IMPLEMENTED (hi);
  return false; // no parent
}

