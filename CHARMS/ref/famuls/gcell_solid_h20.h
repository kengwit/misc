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
#ifndef GCELL_SOLID_H20_H
# define GCELL_SOLID_H20_H

#include "fen.h"
#include "gcell.h"
#include "g3d.h"

class GCELL_SOLID_H20 : public GCELL {

 public: // class members
  
  static const char TYPE_NAME[];

  static const sizet NCHILDREN = 8;

  static const sizet NDETAILFENS = 0; // RAW mid-edge nodes, mid-face nodes, mid-volume node

  static const double _child_map[NCHILDREN][2][3];

  static const sizet NEDGES = 12;

  static const sizet NFACES = 6;
  
 public: // class functions /////////////////////////////////////////

  static bool register_make_func ();

 public: // object functions ////////////////////////////////////////

  /**
   */
  GCELL_SOLID_H20 (vector <FEN *> fens);
  ~GCELL_SOLID_H20 () { }

  /**
     What is the manifold dimension of the gcell?
     (Point=0, line=1, surface=2, or solid=3.) 
   */
  GCELL_MANIFOLD_DIM manifold_dim () { return _conn.manifold_dim (); }

  virtual const string type_name () { return string (TYPE_NAME); }

  /**
     Return connectivity.
  */
  const CONN_BASE *conn () { return &_conn; }

  bool eval_bfun_set (EVALPT *evalpt);

  void map_fen (FEN *fen, POINT *param_loc);
  /**
     Map a parametric point xi, eta, theta (one, two, or all three of these
     coordinates may be used) into a child of this gcell.
  */
  bool map_to_child (POINT &param_loc, GCELL **child, POINT *child_param_loc);
  /**
    Given a parametric location, return the parent and the parametric coordinates
     to which this location maps in the parent of the this cell.
     If the cell has no parent, false is returned and the output parameters are undefined.
     If successful, *parent is the parent of the cell, and
     true is returned.
  */ 
  bool map_to_parent (POINT &param_loc, GCELL **parent, POINT *parent_param_loc);
  /**
     Map fen to a fen in the parent.
  */
  FEN * map_to_parent_fen (FEN *fen);
  /**
     Physically divide the gcell by creating its children (and all the nodes
     these children connect on the higher level).
   */
  void divide (class REF_CTX *ref_ctx);
  /**
     Return the set of detail functions that are needed for the refinement
     of the function associated with the node on input.
  */

  void detail_set (FEN *fen_of_bfun_to_refine, set <FEN *> *rf);
  /**
     Return the set of all refinement functions (i.e. detail functions plus
     the one private function) that are needed for the refinement
     of the function associated with the node on input.
  */

  void complete_refinement_set (FEN *fen_of_bfun_to_refine, set <FEN *> *rf);
  /**
     Return the parent of this gcell.  May be null.
  */
  GCELL *parent () { return _parent;  }
  /**
     Get the vector of the children of the cell.
  */
  GCELL *child (sizet i) const { return _child[i]; }
  /**
     Get the number of children.
  */
  sizet nchildren () const {
    if (_child[0]) return NCHILDREN; // We're not checking every child: if there is one, there have to be all of them
    else           return 0;
  }
  /**
     Had this gcell been divided into children?
  */
  bool divided () const { return (nchildren() > 0); }

  /**
   */
  void debug_display ();
#if defined(GIFACE) && GIFACE
  /**
   */
  GraphicObj *make_elixir_obj () {
    WCRec points[8];
    for (sizet j = 0; j < 8; j++) {
      points[j].x = _conn.fen(j)->ref_loc()(0);
      points[j].y = _conn.fen(j)->ref_loc()(1);
      points[j].z = _conn.fen(j)->ref_loc()(2);
    }
    GraphicObj *go = CreateHexahedron (points);
    return go;
  }

  GraphicObj *make_elixir_obj (vector <FIXED_VECTOR<3> > &p) {
    CHECK(p.size() == 8, EXCEPTION_ILLEGAL_USE, ;);
    WCRec points[8];
    for (sizet j = 0; j < 8; j++) {
      points[j].x = p[j](0);
      points[j].y = p[j](1);
      points[j].z = p[j](2);
    }
    GraphicObj *go = CreateHexahedron (points);
    return go;
  }
#endif

 private: // object data ////////////////////////////////////////////

  CONN_SOLID_20    _conn;
  GCELL_SOLID_H20 *_parent;
  GCELL_SOLID_H20 *_child[NCHILDREN];

 private: // object auxiliary ///////////////////////////////////////

  void eval_ns (GCELL_SOLID_H20 *child, EVALPT *evalpt, double xi, double eta, double theta,
                double dxiparent_dxichild);

};

#endif
