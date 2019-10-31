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
#include "famuls.h"
#include "gcell_line_l2.h"
#include "field.h"
#include "evalpt.h"
#include "ref_ctx.h"

const char GCELL_LINE_L2::TYPE_NAME[] = "line_l2";

const GCELL_LINE_L2::child_param_coord_t GCELL_LINE_L2::_child_map[GCELL_LINE_L2::NCHILDREN] = {
    { { 0}, {-1} },
    { {+1}, { 0} }
  };

GCELL_LINE_L2::GCELL_LINE_L2 (vector <FEN *> fens) : GCELL ()
{
  const sizet nfens = 2;
  CHECK (fens.size () == nfens, EXCEPTION_BAD_VALUE ,;);
  for (unsigned int i = 0; i < nfens; i++) _conn.fen(i) = fens[i];
  _parent = 0;
  for (unsigned int i = 0; i < NCHILDREN; i++) _child[i] = 0;
}

static GCELL *
make (vector <FEN *> fens)
{
  return (new GCELL_LINE_L2 (fens));
}

bool
GCELL_LINE_L2::register_make_func ()
{
  return GCELL::register_make_func (make, string(TYPE_NAME), string ("default"));
}

void
GCELL_LINE_L2::debug_display ()
{
  cerr << "GCELL_LINE_L2 " << this->conn()->fen(0)->id() << " " << this->conn()->fen(1)->id();
}


void
GCELL_LINE_L2::eval_ns (GCELL_LINE_L2 *child, EVALPT *evalpt, double xi, 
                        double dxiparent_dxichild)
{
  const sizet nbfuns = evalpt->nbfuns();
  if (child != 0) {
    for (sizet j = 0; j < GCELL_LINE_L2::NCHILDREN; j++) {
      if (child == this->child(j)) {
        xi    = (1 - xi) / 2 * _child_map[j].lo.xi + (1 + xi) / 2 * _child_map[j].hi.xi;
        //        cerr << "child " << j << " " << " xi " << xi << endl;
        break; // done
      }
    }
    dxiparent_dxichild /= 2;
  } 
  // Evaluate your own basis functions
  // cerr << "xi used " << xi << endl;
  const sizet nfens = 2;
  double N[nfens];
  G3D_line_N (N, xi);
  double N_xi[nfens];
  G3D_line_N_xi (N_xi, xi);
  // Now shove the results into the buffer
  sizet had = 0;
  for (sizet j = 0; j < nbfuns; j++) {
    if ((*evalpt).fen(j) == _conn.fen(0)) {
      evalpt->set_N (j, N[0]);
      evalpt->set_N_der (j, 0, N_xi[0] * dxiparent_dxichild);
      had++;
    }
    if ((*evalpt).fen(j) == _conn.fen(1)) {
      evalpt->set_N (j, N[1]);
      evalpt->set_N_der (j, 0, N_xi[1] * dxiparent_dxichild);
      had++;
    }
    if (had == nfens) break;
  }

  // If you're the child of another cell, ask it to do its job
  if (this->_parent != 0) this->_parent->eval_ns (this, evalpt, xi, dxiparent_dxichild);
}

bool
GCELL_LINE_L2::eval_bfun_set (EVALPT *evalpt)
{
  POINT param_loc = evalpt->param_loc();
  eval_ns (0, evalpt, param_loc(0), 1.0);
  return true;
}

void
GCELL_LINE_L2::divide (REF_CTX *ref_ctx)
{
  if (this->nchildren () == 0) { 
    // Generate refinement nodes first
    vector <FEN *> new_fens(3);
    CONN_POINT_1 vertex_conn[2];
    // Vertex connectivities
    vertex_conn[0].fen(0) = _conn.fen(0);
    vertex_conn[1].fen(0) = _conn.fen(1);
    // Get the refinement nodes for the vertices
    vector <POINT> ref_locs(1);
    ref_locs[0] = _conn.fen(0)->ref_loc();
    new_fens[0] = ref_ctx->get_ref_fen (&vertex_conn[0], 0, ref_locs);
    ref_locs[0] = _conn.fen(1)->ref_loc();
    new_fens[2] = ref_ctx->get_ref_fen (&vertex_conn[1], 0, ref_locs);
    // Get the refinement node for the interior
    POINT midpoint; midpoint.assign (0.5, _conn.fen(0)->ref_loc(), 0.5, _conn.fen(1)->ref_loc());
    ref_locs[0] = midpoint;
    new_fens[1] = ref_ctx->get_ref_fen (&_conn, 0, ref_locs);
    CHECK ((new_fens[0] != 0 &&
            new_fens[1] != 0 &&
            new_fens[2] != 0), EXCEPTION_NULL_PTR ,;);
    // Now generate the new gcells
    vector <FEN *> fens(2);
    fens[0] = new_fens[0];
    fens[1] = new_fens[1];
    {
      GCELL_LINE_L2 *new_gcell = new GCELL_LINE_L2 (fens);
      new_gcell->_parent = this;
      _child[0] = new_gcell;
      this->gcell_group ()->add (new_gcell);
    }
    fens[0] = new_fens[1];
    fens[1] = new_fens[2];
    {
      GCELL_LINE_L2 *new_gcell = new GCELL_LINE_L2 (fens);
      new_gcell->_parent = this;
      _child[1] = new_gcell;
      this->gcell_group ()->add (new_gcell);
    }
  }
}



void
GCELL_LINE_L2::detail_set (FEN *fen_of_bfun_to_refine,set <FEN*> *rf)
{
  int local_index = _conn.local_index (fen_of_bfun_to_refine);
  if (local_index >= 0) {  // Insert nodes
    if (local_index == 1) {
      (*rf).insert (_child[1]->conn()->fen(0));
    } else {
      (*rf).insert (_child[0]->conn()->fen(1));
    }
  }
}



void
GCELL_LINE_L2::complete_refinement_set (FEN *fen_of_bfun_to_refine,set <FEN*> *rf)
{
  int local_index = _conn.local_index (fen_of_bfun_to_refine);
  if (local_index >= 0) {  // Insert nodes
    if (local_index == 1) {
      (*rf).insert (_child[1]->conn()->fen(0));
      (*rf).insert (_child[1]->conn()->fen(1));
    } else {
      (*rf).insert (_child[0]->conn()->fen(0));
      (*rf).insert (_child[0]->conn()->fen(1));
    }
  }
}


/**
   Map the fen to parametric coordinates.  If the field pointer given as argument is
   non-null, try to find the highest gcell that is used by the
   basis function set of that field.
*/
void
GCELL_LINE_L2::map_fen (FEN *fen, POINT *param_loc)
{
  *param_loc = 0;
  if (fen->id() == _conn.fen(0)->id()) {
    (*param_loc)(0) = -1;
    return;
  }
  if (fen->id() == _conn.fen(1)->id()) {
    (*param_loc)(0) = +1;
    return;
  }
  throw EXCEPTION_BAD_VALUE();  // the node is unknown to this cell; raise hell
}

/**
   Given a parametric location, return the child and the parametric coordinates
   to which this location maps in the child of the this cell.
   If the cell has no children, false is returned.  
   If successful, *child is the child into which the node maps as *child_fen, and
   true is returned.
*/
bool
GCELL_LINE_L2::map_to_child (POINT &param_loc, GCELL **child, POINT *child_param_loc)
{
  double xi = param_loc(0);
  if (xi < 0) {
    if (_child[0]) {
      *child = _child[0];
      (*child_param_loc)(0) = 2 * xi + 1; (*child_param_loc)(1) = 0; (*child_param_loc)(2) = 0;
      return true;
    } 
  } else {
    if (_child[1]) {
      *child = _child[1];
      (*child_param_loc)(0) = 2 * xi - 1; (*child_param_loc)(1) = 0; (*child_param_loc)(2) = 0;
      return true;
    } 
  }
  return false; // no children
}


bool
GCELL_LINE_L2::map_to_parent (POINT &param_loc, GCELL **parent, POINT *parent_param_loc)
{
  if (_parent) {
    double xi = param_loc(0);
    for (sizet j = 0; j < GCELL_LINE_L2::NCHILDREN; j++) {
      if (this == _parent->child(j)) {
        (*parent_param_loc)(0) = (1 - xi) / 2 * _child_map[j].lo.xi + (1 + xi) / 2 * _child_map[j].hi.xi;
        (*parent_param_loc)(1) = 0;
        (*parent_param_loc)(2) = 0;
        *parent = _parent;
        return true;
      }
    }
  }
  return false; // no parent
}

