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
#ifndef SELCRIT_FEN_H
# define SELCRIT_FEN_H

#include "selcrit_base.h"

class SELCRIT_FEN : public SELCRIT_BASE {

 public:

  SELCRIT_FEN (sizet id) : _id(id) {}
  ~SELCRIT_FEN () {}

  bool match (FEN *fen) {
    return (fen->id() == _id);
  }

  G3D_box_t box (GMESH *gmesh) { 
    FEN *fen = gmesh->find_fen (_id);
    CHECK (fen, EXCEPTION_BAD_VALUE,;);
    POINT _s = fen->ref_loc ();
    G3D_loc_t p;
    p.x = _s(0); p.y = _s(1); p.z = _s(2);
    G3D_box_t box;
    G3D_INIT_BOX (box, p);
    return box;
  }
  

 private:

  sizet _id;
  
};

#endif
