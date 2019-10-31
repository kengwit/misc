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
#ifndef MAT_MGR_H
# define MAT_MGR_H

#include <list>
#include <iterator>
#include "field.h"
#include "gcell.h"
#include "db.h"
#include "mat.h"
#include "mat_mgr_mat_map.h"

class MAT_MGR {

 private: // class data //////////////////////////////////////////////

  static  int _instance_count;

 public: // object functions ////////////////////////////////////////

  /**
     Construct a mat manager.  Only a single instance may exist
     in the program.  It loads its parameters from the database.
  */
  MAT_MGR (class MGR *mgr);
  /**
     Return the mat named as argument.  The manager consults the database, and
     if the named mat is found, it is loaded and returned; otherwise an exception
     is thrown.
  */
  MAT *mat (string mattype, string matname);

 private: // object data //////////////////////////////////////////

  MGR            *_mgr;
  MAT_MGR_MAT_MAP _mm;
  
};

#endif
