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
#include "bfun_set.h"
#include "load_mgr.h"
#include "mgr.h"

int LOAD_MGR::_instance_count =  0;

LOAD_MGR::LOAD_MGR (MGR *mgr) : _mm (mgr->db())
{
  _mgr = mgr;
  CHECK (_mgr != 0, EXCEPTION_NULL_PTR,;);
  LOAD_MGR::_instance_count++;
  CHECK (LOAD_MGR::_instance_count == 1, EXCEPTION_ILLEGAL_USE,;);
}

LOAD *LOAD_MGR::load (string loadtype, string loadname)
{
  //cerr << "in LOAD_MGR::load: asking for " << loadtype << " " << loadname << endl;
  LOAD *m = _mm.load (loadtype, loadname);
  //  if (m) { cerr << "Produced load " << endl; }
  //  else { cerr << "Sucked" << endl; }
  return m;
}


