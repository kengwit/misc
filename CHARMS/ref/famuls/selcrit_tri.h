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
#ifndef SELCRIT_TRI_H
# define SELCRIT_TRI_H

#include "selcrit_base.h"

class SELCRIT_TRI : public SELCRIT_BASE {

 public:

  SELCRIT_TRI (POINT p, POINT q, POINT r, double dist) : _p(p), _q(q), _r(r), _dist(dist) {}
  ~SELCRIT_TRI () {}

  bool match (FEN *fen) {
    G3D_loc_t p[3];
    p[0].x = _p(0);
    p[0].y = _p(1);
    p[0].z = _p(2);
    p[1].x = _q(0);
    p[1].y = _q(1);
    p[1].z = _q(2);
    p[2].x = _r(0);
    p[2].y = _r(1);
    p[2].z = _r(2);
    G3D_plane_t plane;
    CHECK (G3D_plane_consts_3 (&plane, p+0, p+1, p+2), EXCEPTION_BAD_VALUE,;);
    POINT ref_loc = fen->ref_loc ();
    G3D_loc_t rl; rl.x = ref_loc(0); rl.y = ref_loc(1); rl.z = ref_loc(2);
    double d = fabs (G3D_signed_dist_from_plane (&plane, &rl));
    if (d <= _dist) {
      G3D_tri_t t[1];
      t[0].loc[0] = p[0]; t[0].loc[1] = p[1]; t[0].loc[2] = p[2];
      if (G3D_pt_inside_tri (&rl, &t[0], &plane)) return true;
      {
        POINT r = fen->ref_loc ();
        POINT rq = r - _q;
        POINT rp = r - _p;
        POINT d = _q - _p;
        if (rq.l2_norm() + rp.l2_norm() - d.l2_norm () <= _dist) return true;
      }
      {
        POINT r = fen->ref_loc ();
        POINT rq = r - _r;
        POINT rp = r - _q;
        POINT d = _r - _q;
        if (rq.l2_norm() + rp.l2_norm() - d.l2_norm () <= _dist) return true;
      }
      {
        POINT r = fen->ref_loc ();
        POINT rq = r - _p;
        POINT rp = r - _r;
        POINT d = _p - _r;
        if (rq.l2_norm() + rp.l2_norm() - d.l2_norm () <= _dist) return true;
      }
    }
    return false;
  }

  G3D_box_t box (GMESH *gmesh) { 
    G3D_box_t box;
    G3D_loc_t p;
    p.x = _p(0); p.y = _p(1); p.z = _p(2);
    G3D_INIT_BOX (box, p);
    p.x = _q(0); p.y = _q(1); p.z = _q(2);
    G3D_UPDT_BOX (box, p);
    p.x = _r(0); p.y = _r(1); p.z = _r(2);
    G3D_UPDT_BOX (box, p);
    G3D_INFLATE_BOX (box, _dist);
    return box;
  }

 private:

  POINT _p, _q, _r;
  double _dist;
  
};

#endif
