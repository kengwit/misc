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
#ifndef LOAD_ON_GCELL_H
# define LOAD_ON_GCELL_H

#include "load.h"
#include "point.h"
#include "func.h"

class LOAD_ON_GCELL : public LOAD {
  
 public: // class data //////////////////////////////////////////////

  static const string TYPE_NAME;

 public: // object methods //////////////////////////////////////////

  LOAD_ON_GCELL (DB *db, string name) : LOAD (name) {}
  virtual ~LOAD_ON_GCELL () {}

  string type_name () const { return TYPE_NAME; }

};

#endif

