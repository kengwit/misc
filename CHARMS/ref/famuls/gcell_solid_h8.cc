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
#include "gcell_solid_h8.h"
#include "field.h"
#include "evalpt.h"
#include "ref_ctx.h"

const char GCELL_SOLID_H8::TYPE_NAME[] = "solid_h8";

const double GCELL_SOLID_H8::_child_map[NCHILDREN][2][3] = {
  { {-1,-1,-1}, {0,0,0} },
  { { 0,-1,-1}, {1,0,0} },
  { { 0, 0,-1}, {1,1,0} },
  { {-1, 0,-1}, {0,1,0} },
  { {-1,-1, 0}, {0,0,1} },
  { { 0,-1, 0}, {1,0,1} },
  { { 0, 0, 0}, {1,1,1} },
  { {-1, 0, 0}, {0,1,1} }
};

extern "C" {
#include "g3dhex.h"
}

GCELL_SOLID_H8::GCELL_SOLID_H8 (vector <FEN *> fens) : GCELL ()
{
  sizet nfens = _conn.nfens ();
  CHECK (fens.size () == nfens, EXCEPTION_BAD_VALUE ,;);
  for (unsigned int i = 0; i < nfens; i++) _conn.fen(i) = fens[i];
  _parent = 0;
  for (sizet i = 0; i < NCHILDREN; i++) _child[i] = 0;
}


static GCELL *
make (vector <FEN *> fens)
{
  return (new GCELL_SOLID_H8 (fens));
}

bool
GCELL_SOLID_H8::register_make_func ()
{
  return GCELL::register_make_func (make, string(TYPE_NAME), string ("default"));
}

void
GCELL_SOLID_H8::eval_ns (GCELL_SOLID_H8 *child, EVALPT *evalpt, 
                         double xi, double eta, double theta, double dxiparent_dxichild)
{
  const sizet nbfuns = evalpt->nbfuns();
  if (child != 0) {
    for (sizet j = 0; j < GCELL_SOLID_H8::NCHILDREN; j++) {
      if (child == this->child(j)) {
        const int LO = 0, HI = 1;
        xi    = (1 -    xi) / 2 * _child_map[j][LO][0] + (1 +    xi) / 2 * _child_map[j][HI][0];
        eta   = (1 -   eta) / 2 * _child_map[j][LO][1] + (1 +   eta) / 2 * _child_map[j][HI][1];
        theta = (1 - theta) / 2 * _child_map[j][LO][2] + (1 + theta) / 2 * _child_map[j][HI][2];
        break; // done
      }
    }
    dxiparent_dxichild /= 2;
  } 

  // Evaluate your own basis functions
  const sizet nfens = 8;
  double N[nfens];
  G3D_hex_N (N, xi, eta, theta);
  double N_xi[nfens], N_eta[nfens], N_theta[nfens];
  G3D_hex_N_xi    (   N_xi, xi, eta, theta);
  G3D_hex_N_eta   (  N_eta, xi, eta, theta);
  G3D_hex_N_theta (N_theta, xi, eta, theta);
  // Now shove the results into the buffer
  sizet had = 0;
  for (sizet j = 0; j < nbfuns; j++) {
    for (sizet k = 0; k < nfens; k++) {
      if (evalpt->fen(j) == _conn.fen(k)) {
        evalpt->set_N (j, N[k]);
        evalpt->set_N_der (j, 0, N_xi[k] * dxiparent_dxichild);
        evalpt->set_N_der (j, 1, N_eta[k] * dxiparent_dxichild);
        evalpt->set_N_der (j, 2, N_theta[k] * dxiparent_dxichild);
        had++;
      }
    }
    if (had == nfens) break;
  }

  // If you're the child of another cell, ask it to do its job
  if (this->_parent != 0) this->_parent->eval_ns (this, evalpt, xi, eta, theta, dxiparent_dxichild);
}

bool
GCELL_SOLID_H8::eval_bfun_set (EVALPT *evalpt)
{
  POINT param_loc = evalpt->param_loc();
  eval_ns (0, evalpt, param_loc(0), param_loc(1), param_loc(2), 1.0);
  return true;
}

void
GCELL_SOLID_H8::divide (REF_CTX *ref_ctx)
{
  if (this->nchildren () == 0) { // Not refined yet
    const sizet nfens = 8;
#   define EVN(e,v) g3dhex_vertex_on_edge[e][v]
#   define FVN(f,v) g3dhex_vertex_on_face[f][v]
    // Get the refinement nodes for the vertices
    vector <FEN *> v_fens(nfens);
    CONN_POINT_1 vertex_conn[nfens];
    vector <POINT> ref_locs(1);
    for (sizet i = 0; i < nfens; i++) {
      vertex_conn[i].fen(0) = _conn.fen(i);
      ref_locs[0] = _conn.fen(i)->ref_loc();
      v_fens[i] = ref_ctx->get_ref_fen (&vertex_conn[i], 0, ref_locs);
      CHECK (v_fens[i] != 0, EXCEPTION_NULL_PTR,;);
    }
    // Get the refinement nodes for the edges
    vector <FEN *> e_fens(NEDGES);
    CONN_LINE_2 edge_conn[NEDGES];
    for (sizet i = 0; i < NEDGES; i++) {
      edge_conn[i].fen(0) = _conn.fen(EVN(i,0));
      edge_conn[i].fen(1) = _conn.fen(EVN(i,1));
      FIXED_VECTOR<3> l; l.assign (0.5, _conn.fen(EVN(i,0))->ref_loc(),
                                   0.5, _conn.fen(EVN(i,1))->ref_loc());
      ref_locs[0] = l;
      e_fens[i] = ref_ctx->get_ref_fen (&edge_conn[i], 0, ref_locs);
      CHECK (e_fens[i] != 0, EXCEPTION_NULL_PTR,;);
    }
    // Get the refinement nodes for the faces
    vector <FEN *> f_fens(NFACES);
    CONN_SURF_4 face_conn[NFACES];
    for (sizet i = 0; i < NFACES; i++) {
      FIXED_VECTOR<3> l(0.0);
      for (int j = 0; j < 4; j++) {
        face_conn[i].fen(j) = _conn.fen(FVN(i,j));
        l.add (0.25, _conn.fen(FVN(i,j))->ref_loc());
      }
      ref_locs[0] = l;
      f_fens[i] = ref_ctx->get_ref_fen (&face_conn[i], 0, ref_locs);
      CHECK (f_fens[i] != 0, EXCEPTION_NULL_PTR,;);
    }
    // Get the refinement node for the interior
    POINT midpoint(0.0);
    for (sizet i = 0; i < nfens; i++) {
      midpoint.add (1.0/8, _conn.fen(i)->ref_loc());
    }
    ref_locs[0] = midpoint;
    FEN *m_fen = ref_ctx->get_ref_fen (&_conn, 0, ref_locs);
    CHECK (m_fen != 0, EXCEPTION_NULL_PTR,;);
    // Now generate the new gcells
    vector <FEN *> fens(nfens);
    int nc;
#define V(i) v_fens[i]
#define E(i) e_fens[i]
#define F(i) f_fens[i]
#define M m_fen
#define DEFFEN(n0,n1,n2,n3,n4,n5,n6,n7) {fens[0] = n0; fens[1] = n1; fens[2] = n2; fens[3] = n3; fens[4] = n4; fens[5] = n5; fens[6] = n6; fens[7] = n7;}      
      
    nc = 0; ///////////////////////////////////
    DEFFEN (V(0), E(0), F(0), E(3), E(8), F(1), M, F(4));
    _child[nc] = new GCELL_SOLID_H8 (fens);
    _child[nc]->_parent = this;
    this->gcell_group ()->add (_child[nc]);
    nc = 1; ///////////////////////////////////
    DEFFEN (E(0), V(1), E(1), F(0), F(1), E(9), F(2), M);
    _child[nc] = new GCELL_SOLID_H8 (fens);
    _child[nc]->_parent = this;
    this->gcell_group ()->add (_child[nc]);
    nc = 2; ///////////////////////////////////
    DEFFEN (F(0), E(1), V(2), E(2), M, F(2), E(10), F(3));
    _child[nc] = new GCELL_SOLID_H8 (fens);
    _child[nc]->_parent = this;
    this->gcell_group ()->add (_child[nc]);
    nc = 3; ///////////////////////////////////
    DEFFEN (E(3), F(0), E(2), V(3), F(4), M, F(3), E(11));
    _child[nc] = new GCELL_SOLID_H8 (fens);
    _child[nc]->_parent = this;
    this->gcell_group ()->add (_child[nc]);
    nc = 4; ///////////////////////////////////
    DEFFEN (E(8), F(1), M, F(4), V(4), E(4), F(5), E(7));
    _child[nc] = new GCELL_SOLID_H8 (fens);
    _child[nc]->_parent = this;
    this->gcell_group ()->add (_child[nc]);
    nc = 5; ///////////////////////////////////
    DEFFEN (F(1), E(9), F(2), M, E(4), V(5), E(5), F(5));
    _child[nc] = new GCELL_SOLID_H8 (fens);
    _child[nc]->_parent = this;
    this->gcell_group ()->add (_child[nc]);
    nc = 6; ///////////////////////////////////
    DEFFEN (M, F(2), E(10), F(3), F(5), E(5), V(6), E(6));
    _child[nc] = new GCELL_SOLID_H8 (fens);
    _child[nc]->_parent = this;
    this->gcell_group ()->add (_child[nc]);
    nc = 7; ///////////////////////////////////
    DEFFEN (F(4), M, F(3), E(11), E(7), F(5), E(6), V(7));
    _child[nc] = new GCELL_SOLID_H8 (fens);
    _child[nc]->_parent = this;
    this->gcell_group ()->add (_child[nc]); 
  }
}


void
GCELL_SOLID_H8::detail_set (FEN *fen_of_bfun_to_refine,set <FEN *> *rf)
{
  int local_index = _conn.local_index (fen_of_bfun_to_refine);
  if (local_index >= 0) {  // Insert nodes
    const sizet nfens = 8;
    for (sizet i = 0; i < nfens; i++) {
      if (i != (sizet) local_index) {
        (*rf).insert (_child[local_index]->conn()->fen(i));
      } 
    }
  }
}



void
GCELL_SOLID_H8::complete_refinement_set (FEN *fen_of_bfun_to_refine,set <FEN *> *rf)
{
  int local_index = _conn.local_index (fen_of_bfun_to_refine);
  if (local_index >= 0) {  // Insert nodes
    detail_set (fen_of_bfun_to_refine, rf);
    (*rf).insert(_child[local_index]->conn()->fen(local_index));
  }
}


extern "C" {
#include "g3dhex.h"
}
void
GCELL_SOLID_H8::map_fen (FEN *fen, POINT *param_loc)
{
  for (sizet j = 0; j < _conn.nfens (); j++) {
    if (fen->id() == _conn.fen(j)->id()) {
      (*param_loc)(0) = g3dhex_vertex_param_coord[j][0];
      (*param_loc)(1) = g3dhex_vertex_param_coord[j][1];
      (*param_loc)(2) = g3dhex_vertex_param_coord[j][2];
      return;
    }
  }
  throw EXCEPTION_BAD_VALUE();  // the node is unknown to this cell; raise hell
}

bool
GCELL_SOLID_H8::map_to_child (POINT &param_loc, GCELL **child, POINT *child_param_loc)
{
  if (nchildren() == 0) return false;

  double xi = param_loc(0), eta = param_loc(1), theta = param_loc(2);
  if (xi < 0) {
    if (eta < 0) {
      if (theta < 0) {
        *child = _child[0];
        (*child_param_loc)(0) = 2 * xi    + 1;
        (*child_param_loc)(1) = 2 * eta   + 1;
        (*child_param_loc)(2) = 2 * theta + 1;
      } else {
        *child = _child[4];
        (*child_param_loc)(0) = 2 * xi    + 1;
        (*child_param_loc)(1) = 2 * eta   + 1;
        (*child_param_loc)(2) = 2 * theta - 1;
      }
    } else {
      if (theta < 0) {
        *child = _child[3];
        (*child_param_loc)(0) = 2 * xi    + 1;
        (*child_param_loc)(1) = 2 * eta   - 1;
        (*child_param_loc)(2) = 2 * theta + 1;
      } else {
        *child = _child[7];
        (*child_param_loc)(0) = 2 * xi    + 1;
        (*child_param_loc)(1) = 2 * eta   - 1;
        (*child_param_loc)(2) = 2 * theta - 1;
      }
    }
  } else {
    if (eta < 0) {
      if (theta < 0) {
        *child = _child[1];
        (*child_param_loc)(0) = 2 * xi    - 1;
        (*child_param_loc)(1) = 2 * eta   + 1;
        (*child_param_loc)(2) = 2 * theta + 1;
      } else {
        *child = _child[5];
        (*child_param_loc)(0) = 2 * xi    - 1;
        (*child_param_loc)(1) = 2 * eta   + 1;
        (*child_param_loc)(2) = 2 * theta - 1;
      }
    } else {
      if (theta < 0) {
        *child = _child[2];
        (*child_param_loc)(0) = 2 * xi    - 1;
        (*child_param_loc)(1) = 2 * eta   - 1;
        (*child_param_loc)(2) = 2 * theta + 1;
      } else {
        *child = _child[6];
        (*child_param_loc)(0) = 2 * xi    - 1;
        (*child_param_loc)(1) = 2 * eta   - 1;
        (*child_param_loc)(2) = 2 * theta - 1;
      }
    }
  }
  return true; // done
}


bool
GCELL_SOLID_H8::map_to_parent (POINT &param_loc, GCELL **parent, POINT *parent_param_loc)
{
  if (_parent) {
    double xi = param_loc(0), eta = param_loc(1), theta = param_loc(2);
    for (sizet j = 0; j < GCELL_SOLID_H8::NCHILDREN; j++) {
      if (this == _parent->child(j)) {
        const int LO = 0, HI = 1;
        (*parent_param_loc)(0)
          = (1 -    xi) / 2 * _child_map[j][LO][0] + (1 +    xi) / 2 * _child_map[j][HI][0];
        (*parent_param_loc)(1)
          = (1 -   eta) / 2 * _child_map[j][LO][1] + (1 +   eta) / 2 * _child_map[j][HI][1];
        (*parent_param_loc)(2)
          = (1 - theta) / 2 * _child_map[j][LO][2] + (1 + theta) / 2 * _child_map[j][HI][2];
        *parent = _parent;
        return true;
      }
    }
  }
  return false; // no parent
}

void
GCELL_SOLID_H8::debug_display ()
{
  cerr << "GCELL_SOLID_H8 "
       << this->conn()->fen(0)->id() << " "
       << this->conn()->fen(1)->id() << " "
       << this->conn()->fen(2)->id() << " "
       << this->conn()->fen(3)->id() << " "
       << this->conn()->fen(4)->id() << " "
       << this->conn()->fen(5)->id() << " "
       << this->conn()->fen(6)->id() << " "
       << this->conn()->fen(7)->id();
}

