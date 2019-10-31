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
#ifndef SELCRIT_BOX_H
# define SELCRIT_BOX_H

#include "selcrit_base.h"

class SELCRIT_BOX : public SELCRIT_BASE {

 public:

  SELCRIT_BOX (POINT p, POINT q, double inflate) : _p(p), _q(q), _inflate(inflate) {}
  ~SELCRIT_BOX () {}

  bool match (FEN *fen) {
    G3D_box_t box;
    G3D_loc_t p;
    p.x = _p(0); p.y = _p(1); p.z = _p(2);
    G3D_INIT_BOX (box, p);
    p.x = _q(0); p.y = _q(1); p.z = _q(2);
    G3D_UPDT_BOX (box, p);
    G3D_INFLATE_BOX (box, _inflate);
    G3D_loc_t l;
    POINT ref_loc = fen->ref_loc ();
    l.x = ref_loc(0);
    l.y = ref_loc(1);
    l.z = ref_loc(2);
    return (G3D_PT_INSIDE_BOX (box, l));
  }

  G3D_box_t box (GMESH *gmesh) { 
    G3D_box_t box;
    G3D_loc_t p;
    p.x = _p(0); p.y = _p(1); p.z = _p(2);
    G3D_INIT_BOX (box, p);
    p.x = _q(0); p.y = _q(1); p.z = _q(2);
    G3D_UPDT_BOX (box, p);
    G3D_INFLATE_BOX (box, _inflate);
    return box;
  }

 private:

  POINT _p, _q;
  double _inflate;
  
};

#endif
