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
#include "ecell_xhexvue.h"
#include "evalpt.h"

OFACTORY<ECELL_XHEXVUE, ECELL_XHEXVUE::make_ecell_func, GCELL *> * ECELL_XHEXVUE::_ecell_factory = 0;


bool
ECELL_XHEXVUE::register_make_func (ECELL_XHEXVUE::make_ecell_func make, string gcell_type, string implementation)
{
  if (_ecell_factory == 0) {
    ECELL_XHEXVUE::_ecell_factory = new OFACTORY<ECELL_XHEXVUE, ECELL_XHEXVUE::make_ecell_func, GCELL* >;
  }
  return (ECELL_XHEXVUE::_ecell_factory)->register_make_func (make, gcell_type, implementation);
}

ECELL_XHEXVUE *
ECELL_XHEXVUE::make_ecell (GCELL *gcell, string implementation)
{
  CHECK (ECELL_XHEXVUE::_ecell_factory, EXCEPTION_NULL_PTR,;);
  return (ECELL_XHEXVUE::_ecell_factory)->make (gcell, gcell->type_name (), implementation);
}
