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
#ifndef ALGO_motion_pict_H
# define ALGO_motion_pict_H

#include "fen.h"
#include "algo.h"
#include "gmesh.h"
#include "db.h"
#include "field.h"
#include <stack>

class ALGO_MOTION_PICT : public ALGO {

 public:

  static const int FEN_LAYER = 0;

  static const sizet MAX_COLORS = 13;
  
  static const char *color_names_table[MAX_COLORS]; 

  static const double SHRINK;

 public:

  ALGO_MOTION_PICT(string name, MGR *mgr, GMESH *gmesh);
  ~ALGO_MOTION_PICT ();
  void add_graphics (GCELL *gc, FIELD_VECTOR *geometry);
  void add_graphics (FEN *fen, FIELD_VECTOR *geometry);
  void event_loop (bool stop, string msg = "Stopped: press ctrl-p or push button Proceed");
#if defined (GIFACE) && GIFACE
  void clear() { _gogl.clear(); _gofl.clear(); EMDestroyAllGraphics(ESIModel()); }
#endif
  void set_fen_rep (FEN::FEN_REP fen_rep) {
    _fen_rep = fen_rep;
  }
  void set_level_as_fen_layer (bool b) {
    _level_as_fen_layer = b;
  }
  void set_node_size_mult (double m) {
    _node_size_mult = m;
  }
  double scale () const { return _scale; }
  void set_scale (double scale) { _scale = scale; }
  void update (FIELD_VECTOR *geometry);

 private:

  FEN::FEN_REP _fen_rep;
  double       _node_size_mult;
  double       _node_size;
  bool         _level_as_fen_layer;
  double       _scale;
  
#if defined (GIFACE) && GIFACE
  double node_size (GMESH *gmesh);

 private:

  static bool initialized;

  int gcell_layer;
  int fen_layer;
  EPixel fen_color;
  EPixel gcell_color;
  EPixel edge_color;

  typedef struct {
    GraphicObj *go;
    GCELL *gcell;
  } go_gcell_pair_t;

  list <go_gcell_pair_t > _gogl;

  typedef struct {
    GraphicObj *go;
    FEN *fen;
  } go_fen_pair_t;

  list <go_fen_pair_t > _gofl;
#endif

};



#endif
