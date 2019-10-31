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
#ifndef MGR_H
# define MGR_H

#include <list>
#include <iterator>
#include "field.h"
#include "gcell.h"
#include "db.h"
#include "mesh_mgr.h"
#include "mat_mgr.h"
#include "load_mgr.h"
#include "logger_stream_mgr.h"
#include "logger_stream.h"

class MGR {

 private: // class data //////////////////////////////////////////////

  static  int _instance_count;

 public: // object functions ////////////////////////////////////////

  /**
     Construct a mesh manager.  Only a single instance may exist
     in the program.  It loads its parameters from the database.
  */
  MGR (DB *_db);
  /**
     Destructor.
  */
  ~MGR ();
  /**
     Get the database handle.
  */
  DB *db () const { return _db; }
  /**
     Return the mesh manager component.
  */
  MESH_MGR *mesh_mgr () const { return _mesh_mgr; }
  /**
     Return the material manager component.
  */
  MAT_MGR *mat_mgr () const { return _mat_mgr; }
  /**
     Return the load manager component.
  */
  LOAD_MGR *load_mgr () const { return _load_mgr; }
  /**
     Let the user access the logger stream.
  */
  SMART_HANDLE<LOGGER_STREAM> logger () {
    SMART_HANDLE<LOGGER_STREAM> h (_logger_stream_mgr.logger(""));
    return h;
  }
  /**
     Let the user access the logger stream.
  */
  SMART_HANDLE<LOGGER_STREAM> logger (string where, bool print_exec_time) {
    SMART_HANDLE<LOGGER_STREAM> h (_logger_stream_mgr.logger(where, print_exec_time));
    return h;
  }

 private: // object data //////////////////////////////////////////

  DB                     *_db;
  MESH_MGR               *_mesh_mgr;
  MAT_MGR                *_mat_mgr;
  LOAD_MGR               *_load_mgr;
  LOGGER_STREAM_MGR       _logger_stream_mgr;
  
};

#endif
