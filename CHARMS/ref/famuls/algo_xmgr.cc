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
#include "algo_xmgr.h"
#include "evalpt.h"

ALGO_XMGR::ALGO_XMGR (string name, MGR *mgr)  : ALGO (name, mgr)
{
  _set = 0;
}


void
ALGO_XMGR::run (string legend, FIELD_SCALAR *field, FIELD_VECTOR *geometry)
{
  vector <FIELD_PAIR<1> *> pair_list;
  for (sizet j = 0; j < field->npairs(); j++) {
    FIELD_PAIR<1> *pair = field->field_pair (j);
    pair_list.push_back (pair);
  }
  sort (pair_list.begin (), pair_list.end (), ALGO_XMGR::PAIR_COMPARATOR());

  cout << "# ALGO_XMGR " << endl;
  cout << "#   Field: " << field->name() << endl;
  cout << "@ legend string " << _set << " \"" << legend << "\"" << endl;
  for ( vector <FIELD_PAIR<1> *>::iterator i = pair_list.begin (); i != pair_list.end (); i++) {
    FIXED_VECTOR<1> v = field->evaluate (geometry, (*i)->bfun()->fen());
    FIXED_VECTOR<3> x = geometry->evaluate (geometry, (*i)->bfun()->fen());
    cout << x(0) << " " << v(0) << endl;
  }    
  cout << "&" << endl;
  _set++;
}
