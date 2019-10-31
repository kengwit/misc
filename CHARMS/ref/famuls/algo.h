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
#ifndef ALGO_H
# define ALGO_H

#include <string>
#include "mgr.h"

/**
   This class is of purely organizational nature.
 */
class ALGO {

 public:

  ALGO (std::string name, MGR *mgr) : _name (name), _mgr (mgr) { }
  virtual ~ALGO  () {}
  /**
     Instances of algorithms are generally named.  That should help the user
     to keep track of what is happening where.
  */
  const std::string name () { return _name; }
  /**
     Manager handle.
  */
  MGR *mgr () const { return _mgr; }
  
 private:
  
  std::string _name;
  MGR        *_mgr;
  
};

#endif
