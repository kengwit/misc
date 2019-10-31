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
#include "gmesh.h"
#include <algorithm>

GMESH::GMESH (string name, DB *db)
{
  _name = name;
  _fens.clear ();
  _gsubmeshes.clear  ();
  _fen_map.clear ();
  _max_fen_id = 0;
  G3D_loc_t zero={0,0,0};
  G3D_INIT_BOX (_box, zero);

  // FEN
  string fens_param = db->DB_GET_STRING ("gmeshes/" + _name + "/fen_file");
  CHECK (!fens_param.empty (), EXCEPTION_BAD_VALUE,;);
  CHECK (read_fens (fens_param), EXCEPTION_BAD_VALUE,;);
  // GSUBMESH
  string path = "gmeshes/" + _name + "/gsubmeshes";
  list <string> gsubmeshes_list = db->DB_GET_STRING_LIST (path);
  CHECK (!gsubmeshes_list.empty (), EXCEPTION_BAD_VALUE,;);
  list <string>::const_iterator i = gsubmeshes_list.begin ();
  while (i != gsubmeshes_list.end()) {
    string n = *i;
    GSUBMESH *gsubmesh = new GSUBMESH (n, db, this);
    _gsubmeshes.push_back (gsubmesh);
    i++;
  }
}

extern "C" {
#include "tokensP.h"  
}

bool
GMESH::read_fens (string file)
{
  tokens_parser parser = tokens_new_parser ();
  if (!parser) return false;
  if (!tokens_open_file (parser, (char *) file.c_str ())) return false;
  while (tokens_next_line (parser)) {
    CHECK (tokens_total_of_tokens (parser) == 4, EXCEPTION_BAD_VALUE ,;);
    int id = tokens_token_as_int (parser, 1);
    // The node id is supposed to be positive
    CHECK (id > 0, EXCEPTION_ILLEGAL_USE,;);
    POINT p;
    p[0] = tokens_token_as_double (parser, 2);
    p[1] = tokens_token_as_double (parser, 3);
    p[2] = tokens_token_as_double (parser, 4);
    FEN *fen = new FEN (id, p, this);
    _fen_map.insert (map <int, FEN *>::value_type (fen->id (), fen));
  }
  CHECK (!_fens.empty (), EXCEPTION_BAD_VALUE ,;);
  tokens_close_curr_file (parser);
  tokens_delete_parser(parser);
  return true;
}


void
GMESH::add_fen (FEN *fen)
{
  // Generate an internal id for the node
  // Add to a list
  _fens.push_back (fen);
  // Largest external node id (to be used during refinement)
  _max_fen_id = std::max((unsigned) _max_fen_id, fen->id());
  //cerr << "adding fen " << fen->id() << " at " << fen->ref_loc() << endl;
  G3D_loc_t p;
  p.x = fen->ref_loc()(0);
  p.y = fen->ref_loc()(1);
  p.z = fen->ref_loc()(2);
  if (_fens.size() <= 1) {
    G3D_INIT_BOX (_box, p);
  } else {
    G3D_UPDT_BOX (_box, p);
  }
  //  cerr << "_box " << _box.lower.x << " " << _box.lower.y << " " << _box.lower.z << " " <<
  //    _box.upper.x << " " << _box.upper.y << " " << _box.upper.z << endl;
}

FEN *
GMESH::find_fen (int id)
{
  map <int, FEN *>::iterator i = _fen_map.find (id);
  if (i != _fen_map.end ()) return i->second;
  else                      return 0;
}

GSUBMESH *
GMESH::find_gsubmesh (string name)
{
  gsubmesh_enumerator_t e = gsubmesh_enumerator ();
  GSUBMESH *gsubmesh;
  e.reset ();
  while ((gsubmesh = e.next ())) {
    if (gsubmesh->name () == name)
      return gsubmesh;
  }
  return 0;                     // no such submesh
}

void
GMESH::debug_display ()
{
  cerr << "GMESH " << _name << endl;
  gsubmesh_enumerator_t e = gsubmesh_enumerator ();
  GSUBMESH *gsubmesh;
  e.reset ();
  while ((gsubmesh = e.next ())) {
    gsubmesh->debug_display ();
  }
}

GMESH::~GMESH () {
  // Delete the submeshes
  for (list <GSUBMESH *>::iterator i = _gsubmeshes.begin(); i != _gsubmeshes.end(); i++) {
    delete *i;
  }
  _gsubmeshes.clear();
  // Delete the finite element nodes
  for (list <FEN *>::iterator i = _fens.begin(); i != _fens.end(); i++) {
    delete *i;
  }
  _fens.clear();
  _fen_map.clear();
}
