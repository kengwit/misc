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
#include "fen_maker.h"
#include "gmesh.h"
#include "gsubmesh.h"
#include "ref_ctx.h"
#include "bfun_set.h"
#include "mesh_mgr.h"
#include "fixed_vector.h"
#include "mgr.h"

int MESH_MGR::_instance_count =  0;

MESH_MGR::MESH_MGR (MGR *mgr)
{
  _mgr = mgr;
  CHECK (_mgr != 0, EXCEPTION_NULL_PTR,;);
  MESH_MGR::_instance_count++;
  CHECK (MESH_MGR::_instance_count == 1, EXCEPTION_ILLEGAL_USE,;);
  _mesh_map.clear();
}

GMESH *MESH_MGR::gmesh (string name)
{
  mesh_map::iterator i = _mesh_map.find (name);
  if (i != _mesh_map.end()) {
    return i->second;
  } else {
    GMESH *gmesh = new GMESH (name, _mgr->db());
    CHECK (gmesh, EXCEPTION_NULL_PTR,;);
    _mesh_map.insert (mesh_map::value_type (name, gmesh));
    return gmesh;
  }
}

MESH_MGR::~MESH_MGR () {
  // Delete meshes managed by this guy
  for (mesh_map::iterator i = _mesh_map.begin(); i != _mesh_map.end(); i++) {
    GMESH *gmesh = i->second;
    delete gmesh;
  }
}
