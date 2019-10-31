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
#ifndef FIELD_BASE_H
# define FIELD_BASE_H

#include <algorithm>
#include <list>
#include <vector>
#include <set>
#include <string>
#include "fen.h"
#include "field_pair.h"
#include "bfun_dofparam_pair_id.h"
#include "smart_handle.h"

class BFUN_SET;

class FIELD_BASE {

 public: // declarations ////////////////////////////////////////

  typedef pair<set <BFUN_DOFPARAM_PAIR_ID>::const_iterator,
    set <BFUN_DOFPARAM_PAIR_ID>::const_iterator > PAIRS_ITERATORS;
  
 public: // object functions ////////////////////////////////////////

  FIELD_BASE () : _name(""), _bfun_set(0) {}
  /**
     Make a new field with given name, based on the basis function set.
   */
  FIELD_BASE (string name, SMART_HANDLE<class BFUN_SET> bfun_set) {
    _name = name;
    _bfun_set = bfun_set;
  }

  virtual ~FIELD_BASE () { }

  /**
     Provide the caller with an opaque identifier for a degree of freedom
     associated in the field with given node.  If this identifier is passed
     to *any* field based on the same basis functions set, the field will be
     able to return the value of the degree of freedom in constant time (i.e. quickly). 
   */
  virtual BFUN_DOFPARAM_PAIR_ID dofparam_id (FEN *fen) = 0;
  /**
     Return the basis function set on which this field is based.
  */
  SMART_HANDLE<class BFUN_SET> bfun_set () const { return _bfun_set; }
  /**
     Get the pairs active (ie. non-zero basis function) over a gcell.
  */
  virtual PAIRS_ITERATORS pairs_active_over_gcell (GCELL *gcell) = 0;
  /**
     Is the field active over cell?
  */
  virtual bool active_over (GCELL *) const = 0;
  /**
     Return the bfun associated with the identifier on input.   This operation will be done
     efficiently (in constant time).
  */
  virtual BFUN *bfun (BFUN_DOFPARAM_PAIR_ID dofparam_id) const = 0;
  /**
     Return the name.
  */
  const string name  () const { return _name; }
  
 private: // object data //////////////////////////////////////////

  string                       _name;
  SMART_HANDLE<class BFUN_SET> _bfun_set;

};

#endif
