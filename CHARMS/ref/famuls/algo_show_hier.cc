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
#include "algo_show_hier.h"
#include "args.h"

const double ALGO_SHOW_HIER::SHRINK = 0.9;

#if defined (GIFACE) && GIFACE
const char *ALGO_SHOW_HIER::color_names_table[ALGO_SHOW_HIER::MAX_COLORS] = {
  "black",
  "slateblue",
  "hotpink",
  "lightblue",
  "orange",
  "cyan",
  "blueviolet",
  "yellow",
  "blue",
  "green",
  "blueviolet",
  "red",
  "greenyellow"
}; 

#endif

ALGO_SHOW_HIER::ALGO_SHOW_HIER(string name, MGR *mgr) : ALGO (name, mgr)
{
#if defined (GIFACE) && GIFACE
  ELIXIR_INTERFACE::build_interface ((void *) this, (void (*)(void *, Widget))0);
#endif
  _fen_rep = FEN::MARKER;
  _node_size_mult = 0.06;
  _level_as_fen_layer = true;
}

#if defined (GIFACE) && GIFACE
int
ALGO_SHOW_HIER::fen_layer (BFUN *bfun, map <GCELL *, level_gcell_t> &lgmap)
{
  for (sizet j = 0; j < bfun->ngcells(); j++) {
    GCELL *gcell = bfun->gcell (j);
    map <GCELL *, level_gcell_t>::const_iterator i = lgmap.find (gcell);
    if (i != lgmap.end()) {
      return (i->second)._level;
    }
  }
  return FEN_LAYER;
}

double
ALGO_SHOW_HIER::node_size (GMESH *gmesh)
{
  G3D_box_t mesh_box;
  gmesh->box (&mesh_box);
  double r[3];
  r[0] = G3D_BOX_RANGE (mesh_box, x);
  r[1] = G3D_BOX_RANGE (mesh_box, y);
  r[2] = G3D_BOX_RANGE (mesh_box, z);
  sizet fen_total = gmesh->fen_total();
  double mesh_vol = 1;
  int p = 0;
  for (int i = 0; i < 3; i++) {
    if (r[i] > 0) {
      mesh_vol *= r[i];
      p++;
    }
  }
  double ns = pow (mesh_vol / (fen_total + 1), 1.0/p) * _node_size_mult;
  return ns;
}

#endif // GIFACE
