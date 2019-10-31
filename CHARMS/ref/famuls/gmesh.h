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
#ifndef GMESH_H
# define GMESH_H

#include <list>
#include <vector>
#include <map>
#include "fen.h"
#include "gsubmesh.h"
#include "db.h"
#include "enumerator.h"

class GMESH {

 public: // object functions ////////////////////////////////////////

  typedef ENUMERATOR <list <GSUBMESH *>, GSUBMESH> gsubmesh_enumerator_t;
  typedef ENUMERATOR <list <FEN *>, FEN>           fens_enumerator_t;
  
  /**
   */
  GMESH (string name, DB *db);

  /**
   */
  ~GMESH ();

  /**
     Get the mesh name.
  */
  string name () const { return _name; }
  
  /**
     Get an enumerator of the component submeshes.
  */
  gsubmesh_enumerator_t gsubmesh_enumerator () {
    ENUMERATOR <list <GSUBMESH *>, GSUBMESH> e (_gsubmeshes);
    return e;
  }

  /**
     Find a named submesh.
  */
  GSUBMESH *find_gsubmesh (string gsubmesh_name);
  
  /**
     Get an enumerator of the finite element nodes.
   */
  fens_enumerator_t fens_enumerator () {
    ENUMERATOR <list <FEN *>, FEN> e (_fens);
    return e;
  }
  /**
     Get the total of nodes associated with mesh.
  */
  sizet fen_total () const { return _fens.size(); }
  
  /**
     Find a finite element node by id.
  */
  FEN *find_fen (int id);

  /**
     Add a finite element node.
  */
  void add_fen (FEN *fen);

  /**
   */
  sizet max_fen_id () { return _max_fen_id; }
  /**
   */
  void debug_display ();
  /**
   */
  void box (G3D_box_t *box) const { *box = _box; }
  
 private: // object private methods /////////////////////////////////

  bool read_fens (string file);

 private: // object data ////////////////////////////////////////////

  string                                 _name;
  list <FEN *>                           _fens;
  list <GSUBMESH *>                      _gsubmeshes;
  map <int, FEN *>                       _fen_map;
  sizet                                 _max_fen_id;
  G3D_box_t                              _box;
  
};

#endif
