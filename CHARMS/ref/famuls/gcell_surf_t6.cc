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
#include "gcell_surf_t6.h"
#include "field.h"
#include "evalpt.h"
#include "ref_ctx.h"
extern "C" {
#include "g3dtri6.h"
}

const char GCELL_SURF_T6::TYPE_NAME[] = "surf_t6";

const double GCELL_SURF_T6::refn_param_c[9][2] = {
  {0.25,0}, // 6
  {0.75,0}, // 7
  {0.75,0.25}, // 8
  {0.25,0.75}, // 9
  {0,0.75}, // 10
  {0,0.25}, // 11
  {0.25,0.25}, // 12
  {0.5,0.25}, // 13
  {0.25,0.5} // 14
  };

// The child map records the serial numbers of nodes in
// the refinement of the parent.  Refinement nodes 0..5 are located
// at the nodes of the parent.
const int GCELL_SURF_T6::_child_map[NCHILDREN][6] = {
  {0,3,5,6,12,11},
  {3,1,4,7,8,13},
  {5,4,2,14,9,10},
  {4,5,3,14,12,13}
};

GCELL_SURF_T6::GCELL_SURF_T6 (vector <FEN *> fens) : GCELL () 
{
  const sizet nfens = 6;
  CHECK (fens.size () == nfens, EXCEPTION_BAD_VALUE ,;);
  for (unsigned int i = 0; i < nfens; i++) _conn.fen(i) = fens[i];
  _parent = 0;
  for (unsigned int i = 0; i < NCHILDREN; i++) _child[i] = 0;
}

static GCELL *
make (vector <FEN *> fens)
{
  return (new GCELL_SURF_T6 (fens));
}

bool
GCELL_SURF_T6::register_make_func ()
{
  return GCELL::register_make_func (make, string(TYPE_NAME), string ("default"));
}


void
GCELL_SURF_T6::eval_ns (GCELL_SURF_T6 *child, EVALPT *evalpt,
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
  const sizet nfens = 6;
  double N[nfens];
  G3D_tri6_N (N, xi, eta);
  double N_xi[nfens], N_eta[nfens];
  G3D_tri6_N_r ( N_xi, xi, eta);
  G3D_tri6_N_s (N_eta, xi, eta);
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
GCELL_SURF_T6::eval_bfun_set (EVALPT *evalpt)
{
  POINT param_loc = evalpt->param_loc();
  eval_ns (0, evalpt, param_loc(0), param_loc(1), 1.0);
  return true;
}


void
GCELL_SURF_T6::divide (REF_CTX *ref_ctx)
{
  if (this->nchildren () == 0) { // Not refined yet
    const sizet nfens = 6;
    // Get the refinement nodes for the vertices
    int vix = 0;
    vector <FEN *> v_fens(nfens+9);
    {
      CONN_POINT_1 vertex_conn[nfens];
      vector <POINT> ref_locs(1); 
      for (sizet i = 0; i < nfens; i++) {
        vertex_conn[i].fen(0) = _conn.fen(i);
        ref_locs[0] = _conn.fen(i)->ref_loc();
        v_fens[vix] = ref_ctx->get_ref_fen (&vertex_conn[i], 0, ref_locs);
        CHECK (v_fens[i] != 0, EXCEPTION_NULL_PTR,;);
        vix++;
      }
    }
    // Get the refinement nodes for the outer edges
    {
      CONN_LINE_3 edge_conn;
      vector <POINT> ref_locs(2); 
      for (int i = 0; i < 3; i++) {
        edge_conn.fen(0) = _conn.fen(g3dtri6_vertex_on_edge[i][0]);
        edge_conn.fen(1) = _conn.fen(g3dtri6_vertex_on_edge[i][1]);
        edge_conn.fen(2) = _conn.fen(g3dtri6_vertex_on_edge[i][2]);
        FIXED_VECTOR<3> l(0);
        l.assign (0.5, edge_conn.fen(0)->ref_loc(), 0.5, edge_conn.fen(2)->ref_loc());
        ref_locs[0] = l;
        l.assign (0.5, edge_conn.fen(2)->ref_loc(), 0.5, edge_conn.fen(1)->ref_loc());
        ref_locs[1] = l;
        v_fens[vix] = ref_ctx->get_ref_fen (&edge_conn, 0, ref_locs);
        CHECK (v_fens[vix] != 0, EXCEPTION_NULL_PTR,;);
        vix++;
        v_fens[vix] = ref_ctx->get_ref_fen (&edge_conn, 1, ref_locs);
        CHECK (v_fens[vix] != 0, EXCEPTION_NULL_PTR,;);
        vix++;
      }
    }
    // Now three more refinement nodes in the interior 
    {
      vector <POINT> ref_locs(3); 
      FIXED_VECTOR<3> l(0);
      l.assign (0.5, _conn.fen(5)->ref_loc(), 0.5, _conn.fen(3)->ref_loc()); 
      ref_locs[0] = l;
      l.assign (0.5, _conn.fen(3)->ref_loc(), 0.5, _conn.fen(4)->ref_loc()); 
      ref_locs[1] = l;
      l.assign (0.5, _conn.fen(4)->ref_loc(), 0.5, _conn.fen(5)->ref_loc()); 
      ref_locs[2] = l;
      v_fens[vix] = ref_ctx->get_ref_fen (&_conn, 0, ref_locs);
      CHECK (v_fens[vix] != 0, EXCEPTION_NULL_PTR,;);
      vix++;
      v_fens[vix] = ref_ctx->get_ref_fen (&_conn, 1, ref_locs);
      CHECK (v_fens[vix] != 0, EXCEPTION_NULL_PTR,;);
      vix++;
     v_fens[vix] = ref_ctx->get_ref_fen (&_conn, 2, ref_locs);
      CHECK (v_fens[vix] != 0, EXCEPTION_NULL_PTR,;);
      vix++;
    }
    CHECK (vix == 15, EXCEPTION_BAD_VALUE,;);
    // Now generate the new gcells
    vector <FEN *> fens(6);
    for (int nc = 0; nc < 4; nc++) {
      for (sizet i = 0; i < nfens; i++) {
        fens[i] = v_fens[_child_map[nc][i]];
      }
      _child[nc] = new GCELL_SURF_T6 (fens);
      _child[nc]->_parent = this;
      this->gcell_group ()->add (_child[nc]);
    }
  } 
}



void
GCELL_SURF_T6::complete_refinement_set (FEN *fen_of_bfun_to_refine,set <FEN*> *rf)
{
  int local_index = _conn.local_index (fen_of_bfun_to_refine);
  if (local_index >= 0) {  // Insert nodes
    // The topology is a little bit more complicated than for linear-order elements.
    // Finer function is not activated if its coefficient beta *has* to be equal to zero.
    // That is consistent with the Kronecker-delta algorithm in Krysl, Grinspun, Schroeder 2001.
#undef ACF
# define ACF(whichchild,whichfen) (*rf).insert (_child[whichchild]->conn()->fen(whichfen))

    switch (local_index) {
    case 0: // corner
      ACF(0,0); // private node
      ACF(0,3); ACF(0,5); ACF(1,3); ACF(1,5); ACF(2,3); ACF(2,5);
      break;
    case 1: // corner
      ACF(1,1); // private node
      ACF(1,3); ACF(1,4); ACF(0,3); ACF(0,4); ACF(2,3); ACF(2,4);
      break;
    case 2: // corner
      ACF(2,2); // private node
      ACF(2,4); ACF(2,5); ACF(0,4); ACF(0,5); ACF(1,4); ACF(1,5); 
      break;
    case 3: // mid-edge
      ACF(3,2); // private node
      ACF(0,3); ACF(1,3); ACF(3,3); ACF(3,4); ACF(3,5);
      break;
    case 4: // mid-edge
      ACF(3,0); // private node
      ACF(1,4); ACF(2,4); ACF(3,3); ACF(3,4); ACF(3,5);
      break;
    case 5: // mid-edge
      ACF(3,1); // private node
      ACF(0,5); ACF(2,5); ACF(3,3); ACF(3,4); ACF(3,5);
      break;
    }
  }
//  for (set<FEN*>::iterator i = (*rf).begin(); i != (*rf).end(); i++) { //RAW
//    cout << (*i)->id() <<"\n";  //RAW
//  } //RAW
   
}  



void
GCELL_SURF_T6::detail_set (FEN *fen_of_bfun_to_refine, set <FEN *> *rf)
{
  int local_index = _conn.local_index (fen_of_bfun_to_refine);
  if (local_index >= 0) {  // Insert nodes
    // The topology is a little bit more complicated than for linear-order elements.
    // Finer function is not activated if its coefficient beta *has* to be equal to zero.
    // That is consistent with the Kronecker-delta algorithm in Krysl, Grinspun, Schroeder 2001.
#undef ACF
# define ACF(whichchild,whichfen) (*rf).insert(_child[whichchild]->conn()->fen(whichfen))
    switch (local_index) {
    case 0: // corner
      ACF(0,3); ACF(0,5); ACF(1,3); ACF(1,5); ACF(2,3); ACF(2,5);
      break;
    case 1: // corner
      ACF(1,3); ACF(1,4); ACF(0,3); ACF(0,4); ACF(2,3); ACF(2,4);
      break;
    case 2: // corner
      ACF(2,4); ACF(2,5); ACF(0,4); ACF(0,5); ACF(1,4); ACF(1,5); 
      break;
    case 3: // mid-edge
      ACF(0,3); ACF(1,3); ACF(3,3); ACF(3,4); ACF(3,5);
      break;
    case 4: // mid-edge
      ACF(1,4); ACF(2,4); ACF(3,3); ACF(3,4); ACF(3,5);
      break;
    case 5: // mid-edge
      ACF(0,5); ACF(2,5); ACF(3,3); ACF(3,4); ACF(3,5);
      break;
    }
  }
}  


void
GCELL_SURF_T6::map_fen (FEN *fen, POINT *param_loc)
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
  if (fen->id() == _conn.fen(3)->id()) {
    (*param_loc)(0) = 0.5; (*param_loc)(1) = 0; (*param_loc)(2) = 0;
    return;
  }
  if (fen->id() == _conn.fen(4)->id()) {
    (*param_loc)(0) = 0.5; (*param_loc)(1) = 0.5; (*param_loc)(2) = 0;
    return;
  }
  if (fen->id() == _conn.fen(5)->id()) {
    (*param_loc)(0) = 0; (*param_loc)(1) = 0.5; (*param_loc)(2) = 0;
    return;
  }
  throw EXCEPTION_BAD_VALUE();  // the node is unknown to this cell; raise hell
}

bool
GCELL_SURF_T6::map_to_child (POINT &param_loc, GCELL **child, POINT *child_param_loc)
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
GCELL_SURF_T6::map_to_parent (POINT &param_loc, GCELL **parent, POINT *parent_param_loc)
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

void
GCELL_SURF_T6::debug_display ()
{
  cerr << "GCELL_SURF_T6 " 
       << this->conn()->fen(0)->id() << " "
       << this->conn()->fen(1)->id() << " "
       << this->conn()->fen(2)->id() << " "
       << this->conn()->fen(3)->id() << " "
       << this->conn()->fen(4)->id() << " "
       << this->conn()->fen(5)->id() << " "; 
}

