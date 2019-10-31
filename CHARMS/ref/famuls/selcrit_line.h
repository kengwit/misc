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
#ifndef SELCRIT_LINE_H
# define SELCRIT_LINE_H

#include "selcrit_base.h"

class SELCRIT_LINE : public SELCRIT_BASE {

 public:

  SELCRIT_LINE (POINT p, POINT q, double dist) : _p(p), _q(q), _dist(dist) {}
  ~SELCRIT_LINE () {}

  bool match (FEN *fen) {
    POINT r = fen->ref_loc ();
    POINT rq = r - _q;
    POINT rp = r - _p;
    POINT d = _q - _p;
    return (rq.l2_norm() + rp.l2_norm() - d.l2_norm () <= _dist);
  }

  G3D_box_t box (GMESH *gmesh) { 
    G3D_box_t box;
    G3D_loc_t p;
    p.x = _p(0); p.y = _p(1); p.z = _p(2);
    G3D_INIT_BOX (box, p);
    p.x = _q(0); p.y = _q(1); p.z = _q(2);
    G3D_UPDT_BOX (box, p);
    G3D_INFLATE_BOX (box, _dist);
    return box;
  }

 private:

  POINT _p, _q;
  double _dist;
  
};

#endif
