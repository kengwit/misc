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
#include "load.h"

OFACTORY<LOAD, LOAD::make_load_func, pair <DB *, string> > * LOAD::_load_factory = 0;

bool
LOAD::register_make_func (LOAD::make_load_func make, string loadtype)
{
  if (_load_factory == 0) {
    LOAD::_load_factory = new OFACTORY<LOAD, LOAD::make_load_func, pair <DB*,string> >;
  }
  return (LOAD::_load_factory)->register_make_func (make, loadtype, "default");
}


LOAD *
LOAD::make_load (DB *db, string loadtype, string loadname)
{
  CHECK (LOAD::_load_factory, EXCEPTION_NULL_PTR,;);
  LOAD *m =  (LOAD::_load_factory)->make (pair <DB *, string> (db, loadname), loadtype, "default");
  return m;
}
