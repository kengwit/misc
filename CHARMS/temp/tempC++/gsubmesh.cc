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
#include "gsubmesh.h"
#include "gmesh.h"

GSUBMESH::GSUBMESH (string name, DB *db, GMESH *gmesh)
{
  _name = name;
  _gcell_groups.clear ();
  _gmesh = gmesh;
  // GCELL_GROUP
  string path = "gsubmeshes/" + _name + "/gcell_groups";
  list <string> gcell_group_list = db->DB_GET_STRING_LIST (path);
  CHECK (!gcell_group_list.empty (), EXCEPTION_BAD_VALUE,;);
  list <string>::const_iterator i = gcell_group_list.begin ();
  while (i != gcell_group_list.end()) {
    GCELL_GROUP *gcell_group = new GCELL_GROUP ((*i), db, this);
    _gcell_groups.push_back (gcell_group);
    i++;
  }
}

void
GSUBMESH::add (GCELL_GROUP *gcell_group)
{
  _gcell_groups.push_back (gcell_group);
}

void
GSUBMESH::debug_display ()
{
  cerr << "GSUBMESH " << _name << endl;
  for (list <GCELL_GROUP *>::const_iterator i = _gcell_groups.begin ();
       i != _gcell_groups.end(); i++) {
    GCELL_GROUP *gcell_group = *i;
    gcell_group->debug_display ();
  }
}

GSUBMESH::~GSUBMESH () {
  for (list <GCELL_GROUP *>::const_iterator i = _gcell_groups.begin ();
       i != _gcell_groups.end(); i++) {
    GCELL_GROUP *gcell_group = *i;
    delete gcell_group;
  }
}
