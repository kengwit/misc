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
#ifndef GCELL_LINE_L2_H
# define GCELL_LINE_L2_H

#include "fen.h"
#include "gcell.h"

class GCELL_LINE_L2 : public GCELL {

 public: // class members
  
  static const char TYPE_NAME[];

  static const sizet NCHILDREN = 2;

  static const sizet NDETAILFENS = 1;

  typedef struct param_coord_t {
    double xi;
  }              param_coord_t;
  
  typedef struct child_param_coord_t {
    param_coord_t hi, lo;
  }              child_param_coord_t;

  static const child_param_coord_t _child_map[NCHILDREN];

 public: // class functions /////////////////////////////////////////

  static bool register_make_func ();

 public: // object functions ////////////////////////////////////////

  /**
   */
  GCELL_LINE_L2 (vector <FEN *> fens);
  ~GCELL_LINE_L2 () {}
  /**
     What is the manifold dimension of the gcell?
     (Point=0, line=1, surface=2, or solid=3.) 
   */
  GCELL_MANIFOLD_DIM manifold_dim () { return _conn.manifold_dim (); }

  const string type_name () { return string(TYPE_NAME); }

  /**
     Return connectivity.
  */
  const CONN_BASE *conn () { return &_conn; }
  /**
   */
  bool eval_bfun_set (class EVALPT *evalpt);

  void map_fen (FEN *fen, POINT *param_loc);
    /**
     Physically divide the gcell by creating its children (and all the nodes
     these children connect on the higher level).
   */
  void divide (class REF_CTX *ref_ctx);
  /**
     Return the set of detail functions that are needed for the refinement
     of the function associated with the node on input.
  */
//std::set <FEN *> detail_set (FEN *fen_of_bfun_to_refine);
void detail_set (FEN *fen_of_bfun_to_refine,set <FEN*> *rf );
  /**
     Return the set of all refinement functions (i.e. detail functions plus
     the one private function) that are needed for the refinement
     of the function associated with the node on input.
  */
//  std::set <FEN *>complete_refinement_set (FEN *fen_of_bfun_to_refine);
void complete_refinement_set (FEN *fen_of_bfun_to_refine,set <FEN*> *rf );
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
    if (_child[0] && _child[1]) return NCHILDREN;
    else                        return 0;
  }
  /**
     Had this gcell been divided into children?
  */
  bool divided () const { return (nchildren() > 0); }
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
   */
  void debug_display ();

 private: // object methods /////////////////////////////////////////

  void eval_ns (GCELL_LINE_L2 *child, EVALPT *evalpt, double xi, double dxiparent_dxichild);

 private: // object data ////////////////////////////////////////////

  CONN_LINE_2    _conn;
  GCELL_LINE_L2 *_parent;
  GCELL_LINE_L2 *_child[2];

};

#endif
