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
#ifndef LOAD_MGR_H
# define LOAD_MGR_H

#include <list>
#include <iterator>
#include "field.h"
#include "gcell.h"
#include "db.h"
#include "load.h"
#include "load_mgr_load_map.h"

class LOAD_MGR {

 private: // class data //////////////////////////////////////////////

  static  int _instance_count;

 public: // object functions ////////////////////////////////////////

  /**
     Construct a load manager.  Only a single instance may exist
     in the program.  It loads its parameters from the database.
  */
  LOAD_MGR (class MGR *mgr);
  /**
     Return the load named as argument.  The manager consults the database, and
     if the named load is found, it is loaded and returned; otherwise an exception
     is thrown.
  */
  LOAD *load (string loadtype, string loadname);

 private: // object data //////////////////////////////////////////

  MGR              *_mgr;
  LOAD_MGR_LOAD_MAP _mm;
  
};

#endif
