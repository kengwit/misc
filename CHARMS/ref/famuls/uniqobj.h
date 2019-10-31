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
#ifndef UNIQOBJ_H
# define UNIQOBJ_H

#include "famuls.h"

class UNIQOBJ {

  static sizet _nobj;

 public: // object functions ////////////////////////////////////////

  static sizet nobj () { return _nobj; }

 public: // object functions ////////////////////////////////////////

  UNIQOBJ () { _internal_id = _nobj; _nobj++; }
  ~UNIQOBJ () { }

  sizet id () const { return _internal_id; }
  
 private: // object data ////////////////////////////////////////////
  
  sizet            _internal_id;

};

#endif
