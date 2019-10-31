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
#ifndef SELCRIT_NQUAD_H
# define SELCRIT_NQUAD_H

#include "selcrit_base.h"

class SELCRIT_NQUAD : public SELCRIT_BASE {

 public:

  SELCRIT_NQUAD (int n0, int n1, int n2, int n3, double dist) {
    _n[0] = n0;
    _n[1] = n1;
    _n[2] = n2;
    _n[3] = n3;
    _dist = dist;
  }
  ~SELCRIT_NQUAD () {}

  bool match (FEN *fen) {
    G3D_loc_t p[4];
    GMESH *gmesh = fen->gmesh();
    for (sizet i = 0; i < 4; i++) {
      FEN *afen = gmesh->find_fen (_n[i]);
      CHECK (afen, EXCEPTION_BAD_VALUE,;);
      POINT rl = afen->ref_loc ();
      p[i].x = rl(0);
      p[i].y = rl(1);
      p[i].z = rl(2);
    }
    G3D_plane_t plane;
    CHECK (G3D_plane_consts_4 (&plane, p+0, p+1, p+2, p+3), EXCEPTION_BAD_VALUE,;);
    POINT ref_loc = fen->ref_loc ();
    G3D_loc_t rl; rl.x = ref_loc(0); rl.y = ref_loc(1); rl.z = ref_loc(2);
    double d = fabs (G3D_signed_dist_from_plane (&plane, &rl));
    if (d <= _dist) {
      G3D_tri_t t[2];
      t[0].loc[0] = p[0]; t[0].loc[1] = p[1]; t[0].loc[2] = p[2];
      t[1].loc[0] = p[2]; t[1].loc[1] = p[3]; t[1].loc[2] = p[0]; 
      if (G3D_pt_inside_tri (&rl, &t[0], &plane) ||
          G3D_pt_inside_tri (&rl, &t[1], &plane)) return true;
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
        POINT rq = r - _s;
        POINT rp = r - _r;
        POINT d = _s - _r;
        if (rq.l2_norm() + rp.l2_norm() - d.l2_norm () <= _dist) return true;
      }
      {
        POINT r = fen->ref_loc ();
        POINT rq = r - _p;
        POINT rp = r - _s;
        POINT d = _p - _s;
        if (rq.l2_norm() + rp.l2_norm() - d.l2_norm () <= _dist) return true;
      }
    }
    return false;
  }

  G3D_box_t box (GMESH *gmesh) { 
    if (_n[0] == 0) {
      for (sizet i = 0; i < 4; i++) {
        _n[i] = gmesh->find_fen (_nid[i]);
        CHECK (_n[i], EXCEPTION_BAD_VALUE,;);
        _p[i] = _n[i]->ref_loc ();
      }
    }
    G3D_box_t box;
    G3D_loc_t p;
    p.x = _p(0); p.y = _p(1); p.z = _p(2);
    G3D_INIT_BOX (box, p);
    p.x = _q(0); p.y = _q(1); p.z = _q(2);
    G3D_UPDT_BOX (box, p);
    p.x = _r(0); p.y = _r(1); p.z = _r(2);
    G3D_UPDT_BOX (box, p);
    p.x = _s(0); p.y = _s(1); p.z = _s(2);
    G3D_UPDT_BOX (box, p);
    G3D_INFLATE_BOX (box, _dist);
    return box;
  }

 private:

  POINT _p, _q, _r, _s;
  double _dist;
  
};

#endif
