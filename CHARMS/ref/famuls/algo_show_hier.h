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
#ifndef ALGO_SHOW_HIER_H
# define ALGO_SHOW_HIER_H

#include "fen.h"
#include "algo.h"
#include "gmesh.h"
#include "db.h"
#include "field.h"
#include <stack>

class ALGO_SHOW_HIER : public ALGO {

 public:

  static const int FEN_LAYER = 0;

  static const sizet MAX_COLORS = 13;
  
  static const char *color_names_table[MAX_COLORS]; 

  static const double SHRINK;

 public:

  ALGO_SHOW_HIER(string name, MGR *mgr);
  ~ALGO_SHOW_HIER () {}
  template <int NUM_COMPONENTS>
    void save_as_elixir_file (string file, FIELD<NUM_COMPONENTS> *field);
  void set_fen_rep (FEN::FEN_REP fen_rep) {
    _fen_rep = fen_rep;
  }
  void set_level_as_fen_layer (bool b) {
    _level_as_fen_layer = b;
  }
  void set_node_size_mult (double m) {
    _node_size_mult = m;
  }

 private:

  FEN::FEN_REP _fen_rep;
  double       _node_size_mult;
  bool         _level_as_fen_layer;
  
#if defined (GIFACE) && GIFACE

  class level_gcell_t {
  public:
    level_gcell_t (int level, GCELL *gcell, GraphicObj *gobj) : _level(level), _gcell(gcell), _gobj(gobj) {}
  public:
    int         _level;
    GCELL      *_gcell;
    GraphicObj *_gobj;
  };

  int fen_layer (BFUN *bfun, map <GCELL *, level_gcell_t> &lgmap);
  double node_size (GMESH *gmesh);
#endif

};


template <int NUM_COMPONENTS>
void
ALGO_SHOW_HIER::save_as_elixir_file (string file, FIELD<NUM_COMPONENTS> *field)
{
#if defined (GIFACE) && GIFACE
  FILE *stream = fopen (file.c_str(), "w");
  CHECK (stream, EXCEPTION_NULL_PTR,;);

  GMESH *gmesh = field->bfun_set()->gmesh();

  // Header
  fprintf (stream, "# ELIXIR file (C) 1994 Petr Krysl; version 1\n");

  BOOLEAN suc;
  EPixel fen_color = ColorGetPixelFromString ("red", &suc);
  CHECK (suc, EXCEPTION_NULL_PTR,;);
  EPixel gcell_color = ColorGetPixelFromString ("red", &suc);
  CHECK (suc, EXCEPTION_NULL_PTR,;);
  EPixel edge_color = ColorGetPixelFromString ("white", &suc);
  CHECK (suc, EXCEPTION_NULL_PTR,;);
  // cerr << "White " << fen_color << " Red " << gcell_color << endl;

  map <GCELL *, level_gcell_t> lgmap;

  // Process the gcells
  sizet max_level = 0;
  BFUN_SET::gsubmesh_enumerator_t e = field->bfun_set()->gsubmesh_enumerator ();
  e.reset ();
  GSUBMESH *gsubmesh;
  while (gsubmesh = e.next()) {
    GSUBMESH::gcell_group_enumerator_t gge = gsubmesh->gcell_group_enumerator ();
    gge.reset ();
    GCELL_GROUP *gcell_group;
    while (gcell_group = gge.next()) {
      GCELL_GROUP::gcell_enumerator_t ge = gcell_group->gcell_enumerator ();
      ge.reset ();
      GCELL *gcell;
      while (gcell = ge.next ()) {
        if (gcell->parent () == 0) { // This is the root gcell
          // Now ascend the hierarchy
          stack <level_gcell_t> s;
          s.push (level_gcell_t(1, gcell, 0));
          while (! s.empty ()) {
            level_gcell_t lg = s.top(); s.pop();
            int level = lg._level;
            GCELL *agcell = lg._gcell;
            if (field->active_over (agcell)) {
              GraphicObj *go = agcell->make_elixir_obj ();
              if (go != 0) {
                EASValsSetLayer (level); if (max_level < level) max_level = level;
                EASValsSetFillStyle (FILL_SOLID);
                EASValsSetColor (gcell_color);
                EASValsSetShrink (SHRINK);
                EASValsSetEdgeColor (edge_color);
                EASValsSetEdgeFlag (TRUE);
                EGWithMaskChangeAttributes (ALL_ATTRIB_MASK, go);
              }
              lgmap.insert (map <GCELL *, level_gcell_t>::value_type (agcell, level_gcell_t(level, agcell, go)));
            }
            for (sizet k = 0; k < agcell->nchildren(); k++) s.push (level_gcell_t(level+1, agcell->child(k), 0));
          } // while 
        } // if
      } // while
    } // while
  } // while

  // Re-color gcells: now we know how many levels there are
  for (map <GCELL *, level_gcell_t>::const_iterator i = lgmap.begin(); i != lgmap.end(); i++) {
    level_gcell_t lg = i->second;
    GraphicObj *go = lg._gobj;
    sizet level = lg._level;
    if (go) {
      BOOLEAN suc;
      EPixel color = ColorGetPixelFromString ((char *) color_names_table[level % MAX_COLORS], &suc);
      EASValsSetColor (color);
      EASValsSetEdgeColor (edge_color);
      EASValsSetEdgeFlag (TRUE);
      EGWithMaskChangeAttributes (COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK, go);
      EGGraphicStoreOn (stream, go);
      EGDestroyGraphics (go);
    }
  }
  
  // Finally the nodes
  double node_size = this->node_size (gmesh);
  if (node_size > 0) {
    ESIHandleCmd ("font fixed");
    for (sizet j = 0; j < field->npairs(); j++) {
      BFUN *bfun = field->field_pair(j)->bfun();
      FEN *fen = bfun->fen();
      GraphicObj *go = fen->make_elixir_obj (node_size, _fen_rep);
      EASValsSetMType (FILLED_CIRCLE_MARKER);
      EASValsSetMSize (8);
      if (this->_level_as_fen_layer)
        EASValsSetLayer (this->fen_layer (bfun, lgmap));
      else
        EASValsSetLayer (FEN_LAYER);
      EASValsSetColor (fen_color);
      EASValsSetFillStyle (FILL_SOLID);
      EASValsSetEdgeFlag (0);
      EASValsSetEdgeColor (fen_color);
      EGWithMaskChangeAttributes (MSIZE_MASK | MTYPE_MASK |
                                  COLOR_MASK | LAYER_MASK |
                                  EDGE_FLAG_MASK | FILL_MASK | 
                                  FONT_MASK, go);
      EGGraphicStoreOn (stream, go);
      EGDestroyGraphics (go);
    }
  }

  // Free resources
  fclose (stream);
#endif // GIFACE
}


#endif
