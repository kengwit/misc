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
#ifndef GSUBMESH_H
# define GSUBMESH_H

#include <string>
#include <list>
// #include "gcell.h"
#include "gcell_group.h"
#include "db.h"

class GSUBMESH {

 public: // declarations //////////////////////////////////////////
  
  typedef ENUMERATOR <list <GCELL_GROUP *>, GCELL_GROUP> gcell_group_enumerator_t;

 public: // object functions ////////////////////////////////////////

  /**
   */
  GSUBMESH (string name, DB *db, class GMESH *gmesh);

  /**
     Destructor: the submesh `owns' the gcell groups, and here
     we delete them.
   */
  ~GSUBMESH ();
  /**
   */
  void add (GCELL_GROUP *gcell_group);

  /**
   */
  gcell_group_enumerator_t gcell_group_enumerator () {
    gcell_group_enumerator_t e = gcell_group_enumerator_t (_gcell_groups);
    return e;
  }

  /**
   */
  class GMESH *gmesh () { return _gmesh; }
  
  /**
   */
  const string name () { return _name; }

  /**
   */
  void debug_display ();
  
 private: // object data ////////////////////////////////////////////

  string               _name;
  list <GCELL_GROUP *> _gcell_groups;
  class GMESH         *_gmesh;
  
};

#endif
