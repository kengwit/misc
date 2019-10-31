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
#ifndef SELCRIT_BASE_H
# define SELCRIT_BASE_H

#include "fen.h"
extern "C" {
# include "g3d.h"
}

class SELCRIT_BASE {

 public:

  SELCRIT_BASE () {}
  virtual ~SELCRIT_BASE () {}

  virtual bool match (FEN *fen) = 0;
  virtual G3D_box_t box (GMESH *gmesh) = 0;

 private:

  int _ignore;
  
};

#endif
