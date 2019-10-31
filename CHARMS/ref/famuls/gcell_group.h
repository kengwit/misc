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
#ifndef GCELL_GROUP_H
# define GCELL_GROUP_H

#include <string>
#include <list>
#include "gcell.h"
#include "enumerator.h"
#include "db.h"

class GCELL_GROUP {

 public: // declarations //////////////////////////////////////////
  
  typedef ENUMERATOR <list <GCELL *>, class GCELL> gcell_enumerator_t;

 public: // object functions //////////////////////////////////////

  /**
   */
  GCELL_GROUP (string name, DB *db, class GSUBMESH *gsubmesh);
  /**
     Destructor.  The group `owns' the gcells, and therefore the get deleted
     here.
   */
  ~GCELL_GROUP ();
  /**
   */
  gcell_enumerator_t gcell_enumerator () {
    gcell_enumerator_t e (_gcells);
    return e;
  }
  /**
     Return gsubmesh to which this gcell group belongs.
  */
  class GSUBMESH *gsubmesh () { return _gsubmesh; }
  /**
     
   */
  const string name () { return _name; }
  /**
     Return gcell type as a string.
   */
  const string gcell_type () { return _gcell_type; }

  /**
     Add a gcell to the group.
  */
  void add (GCELL *gcell);
  /**
   */
  void debug_display ();

  
  
 private: // object functions /////////////////////////////////////

  bool read_gcell_file (string file);

 private: // object data //////////////////////////////////////////

  string                _name;
  class GSUBMESH       *_gsubmesh;
  string                _gcell_type;
  list <GCELL *>        _gcells;
  
};

#endif
