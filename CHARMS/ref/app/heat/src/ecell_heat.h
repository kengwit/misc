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
#ifndef ECELL_HEAT_H
# define ECELL_HEAT_H

#include "gcell.h"
#include "ecell.h"

class ECELL_HEAT {

 public:

  ECELL_HEAT () {}
  virtual bool assemble_conductivity_matrix () = 0;
  virtual bool assemble_capacitance_matrix () = 0;

 private:

};

#include "proto_heat.h"

#endif
