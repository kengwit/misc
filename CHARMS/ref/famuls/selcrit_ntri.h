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
#ifndef SELCRIT_NTRI_H
# define SELCRIT_NTRI_H

#include "selcrit_base.h"

class SELCRIT_NTRI : public SELCRIT_BASE {

 public:

  SELCRIT_NTRI (int n0, int n1, int n2, double dist) {
    POINT zero(0);
    _nid[0] = n0; _n[0] = 0; _p[0] = zero;
    _nid[1] = n1; _n[1] = 0; _p[1] = zero;
    _nid[2] = n2; _n[2] = 0; _p[2] = zero;
    _dist = dist;
  }
  ~SELCRIT_NTRI () {}

  bool match (FEN *fen) {
    if (_n[0] == 0) { // Initialize, the first time round
      GMESH *gmesh = fen->gmesh();
      for (sizet i = 0; i < 3; i++) {
        _n[i] = gmesh->find_fen (_nid[i]);
        CHECK (_n[i], EXCEPTION_BAD_VALUE,;);
        _p[i] = _n[i]->ref_loc ();
      }
    }
    // shortcut
    if ((fen == _n[0]) || (fen == _n[1]) || (fen == _n[2])) goto matched;
    // now for the heavy-duty stuff
    {
      G3D_loc_t p[3];
      for (sizet i = 0; i < 3; i++) {
        p[i].x = _p[i](0);
        p[i].y = _p[i](1);
        p[i].z = _p[i](2);
      }
      G3D_plane_t plane;
      CHECK (G3D_plane_consts_3 (&plane, p+0, p+1, p+2), EXCEPTION_BAD_VALUE,;);
      POINT ref_loc = fen->ref_loc ();
      G3D_loc_t rl; rl.x = ref_loc(0); rl.y = ref_loc(1); rl.z = ref_loc(2);
      double d = fabs (G3D_signed_dist_from_plane (&plane, &rl));
      if (d <= _dist) {
        G3D_tri_t t;
        t.loc[0] = p[0]; t.loc[1] = p[1]; t.loc[2] = p[2];
        if (G3D_pt_inside_tri (&rl, &t, &plane)) goto matched;
        {
          POINT rq = ref_loc - _p[0];
          POINT rp = ref_loc - _p[1];
          POINT d = _p[1] - _p[0];
          if (rq.l2_norm() + rp.l2_norm() - d.l2_norm () <= _dist) goto matched;
        }
        {
          POINT rq = ref_loc - _p[1];
          POINT rp = ref_loc - _p[2];
          POINT d = _p[2] - _p[1];
          if (rq.l2_norm() + rp.l2_norm() - d.l2_norm () <= _dist) goto matched;
        }
        {
          POINT rq = ref_loc - _p[2];
          POINT rp = ref_loc - _p[0];
          POINT d = _p[0] - _p[2];
          if (rq.l2_norm() + rp.l2_norm() - d.l2_norm () <= _dist) goto matched;
        }
      }
    }
    return false;
  matched:
    // cerr << "ntri " << _n[0] << " " << _n[1] << " " << _n[2] << " applies to " << fen->id() << endl;
    return true;
  }

  G3D_box_t box (GMESH *gmesh) { 
    if (_n[0] == 0) {
      for (sizet i = 0; i < 3; i++) {
        _n[i] = gmesh->find_fen (_nid[i]);
        CHECK (_n[i], EXCEPTION_BAD_VALUE,;);
        _p[i] = _n[i]->ref_loc ();
      }
    }
    G3D_box_t box;
    G3D_loc_t p;
    p.x = _p[0](0); p.y = _p[0](1); p.z = _p[0](2);
    G3D_INIT_BOX (box, p);
    p.x = _p[1](0); p.y = _p[1](1); p.z = _p[1](2);
    G3D_UPDT_BOX (box, p);
    p.x = _p[2](0); p.y = _p[2](1); p.z = _p[2](2);
    G3D_UPDT_BOX (box, p);
    G3D_INFLATE_BOX (box, _dist);
    return box;
  }


 private:

  int _nid[3];
  FEN *_n[3];
  POINT _p[3];
  double _dist;
  
};

#endif
