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
#ifndef FAMULS_H
# define FAMULS_H

#undef min
#undef max
#include "famexception.h"
#include "db.h"
#include "util.h"
#include <cstdlib>

typedef long unsigned int sizet;

#define NOT_IMPLEMENTED(what) CHECK (!STRINGIZE_(what), EXCEPTION_NOT_IMPLEMENTED,;)


#endif
