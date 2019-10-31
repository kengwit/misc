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
#ifndef GCELL_SOLID_T10_H
# define GCELL_SOLID_T10_H

#include "fen.h"
#include "gcell.h"
#include "g3d.h"

class GCELL_SOLID_T10 : public GCELL {

 public: // class members
  
  static const char TYPE_NAME[];

  static const sizet NCHILDREN = 8;

  static const sizet NDETAILFENS = 25; 

  typedef struct {double B[3][3]; double t[3];} map_struct;

  static const map_struct map_child_to_parent[NCHILDREN];

  static const map_struct map_child_to_parent_reflected[4];

  static const map_struct map_parent_to_child[NCHILDREN];

  static const map_struct map_parent_to_child_reflected[4];

  static const int _child_map[NCHILDREN][10];

  static const int _child_map_reflected[NCHILDREN][10];

  static const double detfen_param_c[25][3];

 public: // class functions /////////////////////////////////////////

  static bool register_make_func ();

 public: // object functions ////////////////////////////////////////

  /**
   */
  GCELL_SOLID_T10 (vector <FEN *> fens);
  ~GCELL_SOLID_T10 () { }

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

  /**
     Evaluate the basis functions.
  */
  bool eval_bfun_set (EVALPT *evalpt);
  /**
     Return the parent of this gcell.  May be null if the
     cell has no parent.
  */
  GCELL *parent () { return _parent; }
  /**
     Get the vector of the children of the cell.
  */
  GCELL *child (sizet i) const { return _child[i]; }
  /**
     Get the number of children.  If the gcell is not refined zero should be
     returned; otherwise the number of child gcells is returned (depending on the type of the gcell).
  */
  sizet nchildren () const {
    if (_child[0]) return NCHILDREN; /* If there is one, there's all of them */
    else           return 0;
  }
  /**
     Had this gcell been divided into children?
  */
  bool divided () const { return (nchildren() > 0); }

  void debug_display ();

  /**
     Physically divide the gcell by creating its children (and all the nodes
     these children connect on the higher level).
   */
  void divide (class REF_CTX *ref_ctx);
  /**
     Return the set of detail functions that are needed for the refinement
     of the function associated with the node on input.
  */
//  std::set <FEN *> detail_set (FEN *fen_of_bfun_to_refine);
void detail_set (FEN *fen_of_bfun_to_refine,set <FEN*> *rf);
  /**
     Return the set of all refinement functions (i.e. detail functions plus
     the one private function) that are needed for the refinement
     of the function associated with the node on input.
  */
//  std::set <FEN *> complete_refinement_set (FEN *fen_of_bfun_to_refine);
  void complete_refinement_set (FEN *fen_of_bfun_to_refine,set <FEN *> *rf );

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
     Refine the gcell by one level of subdivision.  Needs to be defined
     by the implementation of the specialization.
  */
  void refine ();
#if defined(GIFACE) && GIFACE
  /**
   */
  GraphicObj *make_elixir_obj () {
    WCRec points[10];
    for (sizet j = 0; j < _conn.nfens(); j++) {
      points[j].x = _conn.fen(j)->ref_loc()(0);
      points[j].y = _conn.fen(j)->ref_loc()(1);
      points[j].z = _conn.fen(j)->ref_loc()(2);
    }
    GraphicObj *go = CreateTet10 (points);
    return go;
  }
  GraphicObj *make_elixir_obj (vector <FIXED_VECTOR<3> > &p) {
    CHECK(p.size() == _conn.nfens(), EXCEPTION_ILLEGAL_USE,;);
    WCRec points[10];
    for (sizet j = 0; j < _conn.nfens(); j++) {
      points[j].x = p[j](0);
      points[j].y = p[j](1);
      points[j].z = p[j](2);
    }
    GraphicObj *go = CreateTet10 (points);
    return go;
  }

#endif

  bool is_reflection () const { return _is_reflection; }

 private: // object data ////////////////////////////////////////////

  CONN_SOLID_10    _conn;
  GCELL_SOLID_T10 *_parent;
  GCELL_SOLID_T10 *_child[NCHILDREN];
  bool            _is_reflection;

 private: // object auxiliary ///////////////////////////////////////

  void eval_ns (GCELL_SOLID_T10 *child, EVALPT *evalpt,
                double xi, double eta, double theta, double B[3][3]);
  

};

#endif
