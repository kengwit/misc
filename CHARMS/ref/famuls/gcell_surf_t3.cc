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
#include "gcell_surf_t3.h"
#include "field.h"
#include "evalpt.h"
#include "ref_ctx.h"

const char GCELL_SURF_T3::TYPE_NAME[] = "surf_t3";

const int GCELL_SURF_T3::edge_vert_neigh[3][2] = {{0,1},{1,2},{2,0}};
const int GCELL_SURF_T3::vert_edge_neigh[3][2] = {{0,2},{0,1},{1,2}};

GCELL_SURF_T3::GCELL_SURF_T3 (vector <FEN *> fens) : GCELL () 
{
  sizet nfens = _conn.nfens ();
  CHECK (fens.size () == nfens, EXCEPTION_BAD_VALUE ,;);
  for (unsigned int i = 0; i < nfens; i++) _conn.fen(i) = fens[i];
  _parent = 0;
  for (unsigned int i = 0; i < NCHILDREN; i++) _child[i] = 0;
}

static GCELL *
make (vector <FEN *> fens)
{
  return (new GCELL_SURF_T3 (fens));
}

bool
GCELL_SURF_T3::register_make_func ()
{
  return GCELL::register_make_func (make, string(TYPE_NAME), string ("default"));
}


void
GCELL_SURF_T3::eval_ns (GCELL_SURF_T3 *child, EVALPT *evalpt,
                         double xi, double eta, double dxiparent_dxichild)
{
  const sizet nbfuns = evalpt->nbfuns ();
  // RAW: Should be told what needs to be evaluated (connectivity, basis function value, basis
  //      function value + derivatives) in order to be able to optimize the effort
  if (child != 0) {
    if        (child == this->child(0)) {
      xi  = xi / 2;
      eta = eta / 2;
      dxiparent_dxichild /= 2;
    } else if (child == this->child(1)) {
      xi  = xi / 2 + 1.0 / 2;
      eta = eta / 2;
      dxiparent_dxichild /= 2;
    } else if (child == this->child(2)) {
      xi  = xi / 2;
      eta = eta / 2 + 1.0 / 2;
      dxiparent_dxichild /= 2;
    } else if (child == this->child(3)) {
      xi  = -xi / 2 + 1.0 / 2;
      eta = -eta / 2 + 1.0 / 2;
      dxiparent_dxichild /= -2;
    }
  } 

  // Evaluate your own basis functions
  const sizet nfens = 3;
  double N[nfens];
  N[0] = 1 - xi - eta; N[1] = xi; N[2] = eta;
  double N_xi[nfens], N_eta[nfens];
  N_xi[0]  = -1; N_xi[1]  = 1; N_xi[2]  = 0;
  N_eta[0] = -1; N_eta[1] = 0; N_eta[2] = 1;
  // Now shove the results into the buffer
  sizet had = 0;
  for (sizet j = 0; j < nbfuns; j++) {
    for (sizet k = 0; k < nfens; k++) {
      if (evalpt->fen(j) == _conn.fen(k)) {
        evalpt->set_N (j, N[k]);
        evalpt->set_N_der (j, 0, N_xi[k] * dxiparent_dxichild);
        evalpt->set_N_der (j, 1, N_eta[k] * dxiparent_dxichild);
        had++;
      }
    }
    if (had == nfens) break;
  }

  // If you're the child of another cell, ask it to do its job
  if (this->_parent != 0) this->_parent->eval_ns (this, evalpt, xi, eta, dxiparent_dxichild);
}

bool
GCELL_SURF_T3::eval_bfun_set (EVALPT *evalpt)
{
  POINT param_loc = evalpt->param_loc();
  eval_ns (0, evalpt, param_loc(0), param_loc(1), 1.0);
  return true;
}


void
GCELL_SURF_T3::divide (REF_CTX *ref_ctx)
{
  if (this->nchildren () == 0) { // Not refined yet
    const sizet nfens = 3;
    // Get the refinement nodes for the vertices
    vector <FEN *> v_fens(nfens+3);
    CONN_POINT_1 vertex_conn[nfens];
    vector <POINT> ref_locs(1);
    for (sizet i = 0; i < nfens; i++) {
      vertex_conn[i].fen(0) = _conn.fen(i);
      ref_locs[0] = _conn.fen(i)->ref_loc();
      v_fens[i] = ref_ctx->get_ref_fen (&vertex_conn[i], 0, ref_locs);
      CHECK (v_fens[i] != 0, EXCEPTION_NULL_PTR,;);
    }
    // Get the refinement nodes for the edges
    CONN_LINE_2 edge_conn[3];
    for (int i = 0; i < 3; i++) {
      edge_conn[i].fen(0) = _conn.fen(edge_vert_neigh[i][0]);
      edge_conn[i].fen(1) = _conn.fen(edge_vert_neigh[i][1]);
      FIXED_VECTOR<3> l; l.assign (0.5, _conn.fen(edge_vert_neigh[i][0])->ref_loc(),
                                   0.5, _conn.fen(edge_vert_neigh[i][1])->ref_loc());
      ref_locs[0] = l;
      v_fens[i+3] = ref_ctx->get_ref_fen (&edge_conn[i], 0, ref_locs);
      CHECK (v_fens[i+3] != 0, EXCEPTION_NULL_PTR,;);
    }
    // Now generate the new gcells
    vector <FEN *> fens(3);
    int nc;
    nc = 0; ///////////////////////////////////
    fens[0] = v_fens[0];
    fens[1] = v_fens[3];
    fens[2] = v_fens[5];
    _child[nc] = new GCELL_SURF_T3 (fens);
    _child[nc]->_parent = this;
    this->gcell_group ()->add (_child[nc]);
    nc = 1; ///////////////////////////////////
    fens[0] = v_fens[3];
    fens[1] = v_fens[1];
    fens[2] = v_fens[4];
    _child[nc] = new GCELL_SURF_T3 (fens);
    _child[nc]->_parent = this;
    this->gcell_group ()->add (_child[nc]);
    nc = 2; ///////////////////////////////////
    fens[0] = v_fens[5];
    fens[1] = v_fens[4];
    fens[2] = v_fens[2];
    _child[nc] = new GCELL_SURF_T3 (fens);
    _child[nc]->_parent = this;
    this->gcell_group ()->add (_child[nc]);
    nc = 3; ///////////////////////////////////
    fens[0] = v_fens[4];
    fens[1] = v_fens[5];
    fens[2] = v_fens[3];
    _child[nc] = new GCELL_SURF_T3 (fens);
    _child[nc]->_parent = this;
    this->gcell_group ()->add (_child[nc]);
  } 
}



void
GCELL_SURF_T3::detail_set (FEN *fen_of_bfun_to_refine,set <FEN *> *rf)
{
  int local_index = _conn.local_index (fen_of_bfun_to_refine);
  if (local_index >= 0) {
    const sizet nfens = 3;
    for (sizet i = 0; i < nfens; i++) {
      if (i != (sizet) local_index) {
        (*rf).insert (_child[local_index]->conn()->fen(i));
      }
    }
  }
}



void
GCELL_SURF_T3::complete_refinement_set (FEN *fen_of_bfun_to_refine,set <FEN *> *rf)
{
  int local_index = _conn.local_index (fen_of_bfun_to_refine);
  if (local_index >= 0) {
    this->detail_set (fen_of_bfun_to_refine,rf); // details
    (*rf).insert (_child[local_index]->conn()->fen(local_index)); // private fun
  }
}


void
GCELL_SURF_T3::map_fen (FEN *fen, POINT *param_loc)
{
  if (fen->id() == _conn.fen(0)->id()) {
    (*param_loc)(0) = 0; (*param_loc)(1) = 0; (*param_loc)(2) = 0;
    return;
  }
  if (fen->id() == _conn.fen(1)->id()) {
    (*param_loc)(0) = +1; (*param_loc)(1) = 0; (*param_loc)(2) = 0;
    return;
  }
  if (fen->id() == _conn.fen(2)->id()) {
    (*param_loc)(0) = 0; (*param_loc)(1) = +1; (*param_loc)(2) = 0;
    return;
  }
  throw EXCEPTION_BAD_VALUE();  // the node is unknown to this cell; raise hell
}

bool
GCELL_SURF_T3::map_to_child (POINT &param_loc, GCELL **child, POINT *child_param_loc)
{
  if (nchildren() == 0) return false;
  double xi = param_loc(0), eta = param_loc(1);
  if (xi < 0.5) {
    if (eta > 0.5) {
      *child = _child[2];
      (*child_param_loc)(0) = 2 * xi; (*child_param_loc)(1) = 2 * eta - 1; (*child_param_loc)(2) = 0;
    } else {
      if (0.5 - xi - eta > 0) {
        *child = _child[0];
        (*child_param_loc)(0) = 2 * xi; (*child_param_loc)(1) = 2 * eta; (*child_param_loc)(2) = 0;
      } else {
        *child = _child[3];
        (*child_param_loc)(0) = - 2 * xi + 1; (*child_param_loc)(1) = - 2 * eta + 1; (*child_param_loc)(2) = 0;
      }
    }
  } else {
    *child = _child[1];
    (*child_param_loc)(0) = 2 * xi - 1; (*child_param_loc)(1) = 2 * eta; (*child_param_loc)(2) = 0;
  }
  return true; // no children
}


bool
GCELL_SURF_T3::map_to_parent (POINT &param_loc, GCELL **parent, POINT *parent_param_loc)
{
  if (_parent) {
    double xi = param_loc(0), eta = param_loc(1);
    if        (this == _parent->child(0)) {
      (*parent_param_loc)(0) = xi / 2;
      (*parent_param_loc)(1) = eta / 2;
      (*parent_param_loc)(2) = 0;
      *parent = _parent;
      return true;
    } else if (this == _parent->child(1)) {
      (*parent_param_loc)(0) = xi / 2 + 1.0 / 2;
      (*parent_param_loc)(1) = eta / 2;
      (*parent_param_loc)(2) = 0;
      *parent = _parent;
      return true;
    } else if (this == _parent->child(2)) {
      (*parent_param_loc)(0) = xi / 2;
      (*parent_param_loc)(1) = eta / 2 + 1.0 / 2;
      (*parent_param_loc)(2) = 0;
      *parent = _parent;
      return true;
    } else if (this == _parent->child(3)) {
      (*parent_param_loc)(0) = -xi / 2 + 1.0 / 2;
      (*parent_param_loc)(1) = -eta / 2 + 1.0 / 2;
      (*parent_param_loc)(2) = 0;
      *parent = _parent;
      return true;
    }
  }
  return false; // no parent
}
