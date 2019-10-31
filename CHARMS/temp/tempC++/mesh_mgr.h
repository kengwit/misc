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
#ifndef MESH_MGR_H
# define MESH_MGR_H

#include <list>
#include <iterator>
#include "field.h"
#include "gcell.h"
#include "db.h"

class MESH_MGR {

 private: // class data //////////////////////////////////////////////

  static  int _instance_count;

 public: // object functions ////////////////////////////////////////

  /**
     Construct a mesh manager.  Only a single instance may exist
     in the program.  It loads its parameters from the database.
  */
  MESH_MGR (class MGR *mgr);
  /**
     Destructor: delete all managed meshes.
  */
  ~MESH_MGR ();
  /**
     Return the mesh named as argument.  The manager consults the database, and
     if the named mesh is found, it is loaded and returned; otherwise an exception
     is thrown.
  */
  GMESH *gmesh (string name);

 private: // object data //////////////////////////////////////////

  typedef map <std::string, GMESH * > mesh_map;
  class MGR                  *_mgr;
  mesh_map                    _mesh_map;
  
};

#endif
