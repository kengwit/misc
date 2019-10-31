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
#include "mat.h"
#include "bfun_set.h"
#include "mat_mgr.h"
#include "mgr.h"

int MAT_MGR::_instance_count =  0;

MAT_MGR::MAT_MGR (MGR *mgr) : _mm (mgr->db())
{
  _mgr = mgr;
  CHECK (_mgr != 0, EXCEPTION_NULL_PTR,;);
  MAT_MGR::_instance_count++;
  CHECK (MAT_MGR::_instance_count == 1, EXCEPTION_ILLEGAL_USE,;);
}

MAT *MAT_MGR::mat (string mattype, string matname)
{
  //cerr << "in MAT_MGR::mat: asking for " << mattype << " " << matname << endl;
  MAT *m = _mm.mat (mattype, matname);
  //  if (m) { cerr << "Produced material" << endl; }
  //  else { cerr << "Sucked" << endl; }
  return m;
}


