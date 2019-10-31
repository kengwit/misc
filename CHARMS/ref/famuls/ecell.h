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
#ifndef ECELL_H
# define ECELL_H

#include "gcell.h"
#include "watchpoint.h"

class ECELL {

 public:

  ECELL (GCELL *gcell) { _gcell = gcell; }
  virtual ~ECELL () {}
  /**
     Get the associated gcell.
  */
  GCELL *gcell () const { return _gcell; }
  /**
   */
  class ECELL_GROUP *ecell_group () const { return _ecell_group; }
  /**
     Set the ecell group to which this ecell belongs.
  */
  void set_ecell_group (class ECELL_GROUP *ecell_group) { _ecell_group = ecell_group; }

  virtual void set_watch_point_value (WATCH_POINT *watch_point, POINT& param_loc) {/* DO NOTHING */};

 private:

  GCELL *_gcell;
  class ECELL_GROUP *_ecell_group;
  
};

#endif
