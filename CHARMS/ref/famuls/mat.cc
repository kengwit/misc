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

OFACTORY<MAT, MAT::make_mat_func, pair <DB *, string> > * MAT::_mat_factory = 0;

bool
MAT::register_make_func (MAT::make_mat_func make, string mattype)
{
  if (_mat_factory == 0) {
    MAT::_mat_factory = new OFACTORY<MAT, MAT::make_mat_func, pair <DB*,string> >;
  }
  return (MAT::_mat_factory)->register_make_func (make, mattype, "default");
}


MAT *
MAT::make_mat (DB *db, string mattype, string matname)
{
  CHECK (MAT::_mat_factory, EXCEPTION_NULL_PTR,;);
  MAT *m =  (MAT::_mat_factory)->make (pair <DB *, string> (db, matname), mattype, "default");
  return m;
}
