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
#include "gcell_solid_t4.h"
#include "field.h"
#include "evalpt.h"
#include "ref_ctx.h"

const char GCELL_SOLID_T4::TYPE_NAME[] = "solid_t4";

const GCELL_SOLID_T4::map_struct GCELL_SOLID_T4::map_child_to_parent[NCHILDREN] = {
{{{0.5,0,0},{0,0.5,0},{0,0,0.5}},{0,0,0}}, // child 0
{{{0.5,0,0},{0,0.5,0},{0,0,0.5}},{0.5,0,0}}, // child 1
{{{0.5,0,0},{0,0.5,0},{0,0,0.5}},{0,0.5,0}}, // child 2
{{{0.5,0,0},{0,0.5,0},{0,0,0.5}},{0,0,0.5}}, // child 3
{{{-0.5,-0.5,0},{0,0.5,0},{0,0,-0.5}},{0.5,0,0.5}}, // child 4
{{{0,0,-0.5},{0.5,0,0.5},{-0.5,-0.5,0}},{0.5,0,0.5}}, // child 5
{{{0.5,0.5,0},{-0.5,0,0},{0,0,0.5}},{0,0.5,0}}, // child 6
{{{0,0.5,0},{0,-0.5,-0.5},{-0.5,-0.5,0}},{0,0.5,0.5}}, // child 7
};

const GCELL_SOLID_T4::map_struct GCELL_SOLID_T4::map_child_to_parent_reflected[4] = {
{{{-0.5,0,-0.5},{0,0,0.5},{0.5,0.5,0}},{0.5,0,0}}, // child 4
{{{0,-0.5,0},{0.5,0.5,0},{0,0,0.5}},{0.5,0,0}}, // child 5
{{{0,0.5,0.5},{0,0,-0.5},{-0.5,-0.5,0}},{0,0.5,0.5}}, // child 6
{{{0.5,0,0},{-0.5,-0.5,0},{0,0,-0.5}},{0,0.5,0.5}}, // child 7
};

const GCELL_SOLID_T4::map_struct GCELL_SOLID_T4::map_parent_to_child[NCHILDREN] = {
{{{2,-0,-0},{0,2,-0},{0,0,2}},{-0,-0,-0}}, // child 0
{{{2,-0,-0},{0,2,-0},{0,0,2}},{-1,-0,-0}}, // child 1
{{{2,-0,-0},{0,2,-0},{0,0,2}},{-0,-1,-0}}, // child 2
{{{2,-0,-0},{0,2,-0},{0,0,2}},{-0,-0,-1}}, // child 3
{{{-2,-2,-0},{0,2,0},{0,0,-2}},{1,-0,1}}, // child 4
{{{2,2,0},{-2,-2,-2},{-2,0,0}},{-1,2,1}}, // child 5
{{{0,-2,-0},{2,2,-0},{0,0,2}},{1,-1,-0}}, // child 6
{{{-2,0,-2},{2,0,0},{-2,-2,0}},{1,-0,1}}, // child 7
};

const GCELL_SOLID_T4::map_struct GCELL_SOLID_T4::map_parent_to_child_reflected[4] = {
{{{-2,-2,0},{2,2,2},{0,2,0}},{1,-1,-0}}, // child 4
{{{2,2,-0},{-2,0,0},{0,0,2}},{-1,1,-0}}, // child 5
{{{-2,-2,-2},{2,2,0},{0,-2,0}},{2,-1,1}}, // child 6
{{{2,0,0},{-2,-2,-0},{0,0,-2}},{-0,1,1}}, // child 7
};

#define MCP map_child_to_parent
#define MCPR map_child_to_parent_reflected
#define MPC map_parent_to_child
#define MPCR map_parent_to_child_reflected

const int GCELL_SOLID_T4::_child_map[NCHILDREN][4] = {
  {0,4,6,7}, // child 0
  {4,1,5,8}, // child 1
  {6,5,2,9}, // child 2
  {7,8,9,3}, // child 3
  {8,7,9,4}, // child 4
  {8,5,4,9}, // child 5
  {6,4,5,9}, // child 6
  {9,6,4,7}  // child 7
};

const int GCELL_SOLID_T4::_child_map_reflected[NCHILDREN][4] = {
  {0,4,6,7}, // child 0
  {4,1,5,8}, // child 1
  {6,5,2,9}, // child 2
  {7,8,9,3}, // child 3
  {4,7,8,6}, // child 4
  {4,5,6,8}, // child 5
  {9,6,5,8}, // child 6
  {9,8,7,6}  // child 7
};

#undef MAT_3X3
#define MAT_3X3(A_, op, B_)                     \
 {                                              \
   A_[0][0] op B_[0][0];                        \
   A_[0][1] op B_[0][1];                        \
   A_[0][2] op B_[0][2];                        \
   A_[1][0] op B_[1][0];                        \
   A_[1][1] op B_[1][1];                        \
   A_[1][2] op B_[1][2];                        \
   A_[2][0] op B_[2][0];                        \
   A_[2][1] op B_[2][1];                        \
   A_[2][2] op B_[2][2];                        \
 }
#undef MULT_AB_3X3
#define MULT_AB_3X3(C_, A_, B_)                                             \
 {                                                                          \
   C_[0][0]  = (A_[0][0]*B_[0][0] + A_[0][1]*B_[1][0] + A_[0][2]*B_[2][0]); \
   C_[0][1]  = (A_[0][0]*B_[0][1] + A_[0][1]*B_[1][1] + A_[0][2]*B_[2][1]); \
   C_[0][2]  = (A_[0][0]*B_[0][2] + A_[0][1]*B_[1][2] + A_[0][2]*B_[2][2]); \
   C_[1][0]  = (A_[1][0]*B_[0][0] + A_[1][1]*B_[1][0] + A_[1][2]*B_[2][0]); \
   C_[1][1]  = (A_[1][0]*B_[0][1] + A_[1][1]*B_[1][1] + A_[1][2]*B_[2][1]); \
   C_[1][2]  = (A_[1][0]*B_[0][2] + A_[1][1]*B_[1][2] + A_[1][2]*B_[2][2]); \
   C_[2][0]  = (A_[2][0]*B_[0][0] + A_[2][1]*B_[1][0] + A_[2][2]*B_[2][0]); \
   C_[2][1]  = (A_[2][0]*B_[0][1] + A_[2][1]*B_[1][1] + A_[2][2]*B_[2][1]); \
   C_[2][2]  = (A_[2][0]*B_[0][2] + A_[2][1]*B_[1][2] + A_[2][2]*B_[2][2]); \
 }


GCELL_SOLID_T4::GCELL_SOLID_T4 (vector <FEN *> fens) : GCELL ()
{
  sizet nfens = _conn.nfens ();
  CHECK (fens.size () == nfens, EXCEPTION_BAD_VALUE ,;);
  for (unsigned int i = 0; i < nfens; i++) _conn.fen(i) = fens[i];
  _parent = 0;
  for (unsigned int i = 0; i < NCHILDREN; i++) _child[i] = 0;
  _is_reflection = false;
}


static GCELL *
make (vector <FEN *> fens)
{
  return (new GCELL_SOLID_T4 (fens));
}

bool
GCELL_SOLID_T4::register_make_func ()
{
  return GCELL::register_make_func (make, string(TYPE_NAME), string ("default"));
}

void
GCELL_SOLID_T4::eval_ns (GCELL_SOLID_T4 *child, EVALPT *evalpt,
                         double xi, double eta, double theta, double B[3][3])
{
  const sizet nbfuns = evalpt->nbfuns ();
  // RAW: Should be told what needs to be evaluated (connectivity, basis function value, basis
  //      function value + derivatives) in order to be able to optimize the effort
  if (child != 0) {
    for (sizet j = 0; j < GCELL_SOLID_T4::NCHILDREN; j++) {
      if (child == this->child(j)) {
        if ((j < 4) || (! this->_is_reflection)) {
          double newxi    = MCP[j].B[0][0] * xi + MCP[j].B[0][1] * eta + MCP[j].B[0][2] * theta + MCP[j].t[0];
          double neweta   = MCP[j].B[1][0] * xi + MCP[j].B[1][1] * eta + MCP[j].B[1][2] * theta + MCP[j].t[1];
          double newtheta = MCP[j].B[2][0] * xi + MCP[j].B[2][1] * eta + MCP[j].B[2][2] * theta + MCP[j].t[2];
          xi = newxi; eta = neweta; theta = newtheta;
          double prevB[3][3]; MAT_3X3 (prevB, =, B);
          MULT_AB_3X3 (B, MCP[j].B, prevB);
        } else {
          j -= 4;
          double newxi    = MCPR[j].B[0][0] * xi + MCPR[j].B[0][1] * eta + MCPR[j].B[0][2] * theta + MCPR[j].t[0];
          double neweta   = MCPR[j].B[1][0] * xi + MCPR[j].B[1][1] * eta + MCPR[j].B[1][2] * theta + MCPR[j].t[1];
          double newtheta = MCPR[j].B[2][0] * xi + MCPR[j].B[2][1] * eta + MCPR[j].B[2][2] * theta + MCPR[j].t[2];
          xi = newxi; eta = neweta; theta = newtheta;
          double prevB[3][3]; MAT_3X3 (prevB, =, B);
          MULT_AB_3X3 (B, MCPR[j].B, prevB);
        }
        break; // done
      } // if for child
    } // loop thru children
  } 

  // Evaluate your own basis functions
  const sizet nfens = 4;
  double N[nfens]; // BELOW (the basis functions)
  N[0] = 1 - xi - eta - theta; N[1] = xi; N[2] = eta; N[3] = theta; 
  double N_xi[nfens], N_eta[nfens], N_theta[nfens];
  N_xi[0]    = -1; N_xi[1]    = 1; N_xi[2]    = 0; N_xi[3]    = 0;
  N_eta[0]   = -1; N_eta[1]   = 0; N_eta[2]   = 1; N_eta[3]   = 0;
  N_theta[0] = -1; N_theta[1] = 0; N_theta[2] = 0; N_theta[3] = 1;
  // Now shove the results into the buffer
  sizet had = 0;
  for (sizet j = 0; j < nbfuns; j++) {
    for (sizet k = 0; k < nfens; k++) {
      if (evalpt->fen(j) == _conn.fen(k)) {
        evalpt->set_N (j, N[k]);
        for (sizet m = 0; m < 3; m++) {
          evalpt->set_N_der (j, m, B[0][m] * N_xi[k] + B[1][m] * N_eta[k] + B[2][m] * N_theta[k]);
        }
        had++;
      }
    }
    if (had == nfens) break;
  }

  // If you're the child of another cell, ask it to do its job
  if (this->_parent != 0) this->_parent->eval_ns (this, evalpt, xi, eta, theta, B);
}

bool
GCELL_SOLID_T4::eval_bfun_set (EVALPT *evalpt)
{
  double B[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
  POINT param_loc = evalpt->param_loc();
  eval_ns (0, evalpt, param_loc(0), param_loc(1), param_loc(2), B);
  return true;
}

void
GCELL_SOLID_T4::divide (REF_CTX *ref_ctx)
{
  if (this->nchildren () == 0) { // Not refined yet
    const sizet nfens = 4;
#   define EVN(e,v) g3dtet_vertex_on_edge[e][v]
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
    vector <FEN *> e_fens(6);
    CONN_LINE_2 edge_conn[6];
    for (int i = 0; i < 6; i++) {
      edge_conn[i].fen(0) = _conn.fen(EVN(i,0));
      edge_conn[i].fen(1) = _conn.fen(EVN(i,1));
      FIXED_VECTOR<3> l; l.assign (0.5, _conn.fen(EVN(i,0))->ref_loc(),
                                   0.5, _conn.fen(EVN(i,1))->ref_loc());
      ref_locs[0] = l;
      e_fens[i] = ref_ctx->get_ref_fen (&edge_conn[i], 0, ref_locs);
      CHECK (e_fens[i] != 0, EXCEPTION_NULL_PTR,;);
    }
    // Now generate the new gcells
    vector <FEN *> fens(4);
    int nc;
    nc = 0; ///////////////////////////////////
    fens[0] = v_fens[0];
    fens[1] = e_fens[0];
    fens[2] = e_fens[2];
    fens[3] = e_fens[3];
    _child[nc] = new GCELL_SOLID_T4 (fens);
    _child[nc]->_parent = this;
    this->gcell_group ()->add (_child[nc]);
    _child[nc]->_is_reflection = this->_is_reflection;
    nc = 1; ///////////////////////////////////
    fens[0] = e_fens[0];
    fens[1] = v_fens[1];
    fens[2] = e_fens[1];
    fens[3] = e_fens[4];
    _child[nc] = new GCELL_SOLID_T4 (fens);
    _child[nc]->_parent = this;
    this->gcell_group ()->add (_child[nc]);
    _child[nc]->_is_reflection = this->_is_reflection;
    nc = 2; ///////////////////////////////////
    fens[0] = e_fens[2];
    fens[1] = e_fens[1];
    fens[2] = v_fens[2];
    fens[3] = e_fens[5];
    _child[nc] = new GCELL_SOLID_T4 (fens);
    _child[nc]->_parent = this;
    this->gcell_group ()->add (_child[nc]);
    _child[nc]->_is_reflection = this->_is_reflection;
    nc = 3; ///////////////////////////////////
    fens[0] = e_fens[3];
    fens[1] = e_fens[4];
    fens[2] = e_fens[5];
    fens[3] = v_fens[3];
    _child[nc] = new GCELL_SOLID_T4 (fens);
    _child[nc]->_parent = this;
    this->gcell_group ()->add (_child[nc]);
    _child[nc]->_is_reflection = this->_is_reflection;

    // The mapping of the four octehedral children depends on if the parent is a reflection or not.
    if (this-> _is_reflection == false) {
      nc = 4; ///////////////////////////////////
      fens[0] = e_fens[4];
      fens[1] = e_fens[3];
      fens[2] = e_fens[5];
      fens[3] = e_fens[0];
      _child[nc] = new GCELL_SOLID_T4 (fens);
      _child[nc]->_parent = this;
      this->gcell_group ()->add (_child[nc]);
      _child[nc]->_is_reflection = false;
      nc = 5; ///////////////////////////////////
      fens[0] = e_fens[4];
      fens[1] = e_fens[1];
       fens[2] = e_fens[0];
      fens[3] = e_fens[5];
      _child[nc] = new GCELL_SOLID_T4 (fens);
      _child[nc]->_parent = this;
      this->gcell_group ()->add (_child[nc]);
      _child[nc]->_is_reflection = false;
      nc = 6; ///////////////////////////////////
      fens[0] = e_fens[2];
      fens[1] = e_fens[0];
      fens[2] = e_fens[1];
      fens[3] = e_fens[5];
      _child[nc] = new GCELL_SOLID_T4 (fens);
      _child[nc]->_parent = this;
      this->gcell_group ()->add (_child[nc]);
      _child[nc]->_is_reflection = true;
      nc = 7; ///////////////////////////////////
      fens[0] = e_fens[5];
      fens[1] = e_fens[2];
      fens[2] = e_fens[0];
      fens[3] = e_fens[3];
      _child[nc] = new GCELL_SOLID_T4 (fens);
      _child[nc]->_parent = this;
      this->gcell_group ()->add (_child[nc]);
      _child[nc]->_is_reflection = true;
    } else {
      nc = 4; ///////////////////////////////////
      fens[0] = e_fens[0];
      fens[1] = e_fens[3];
      fens[2] = e_fens[4];
      fens[3] = e_fens[2];
      _child[nc] = new GCELL_SOLID_T4 (fens);
      _child[nc]->_parent = this;
      this->gcell_group ()->add (_child[nc]);
      _child[nc]->_is_reflection = false;
      nc = 5; ///////////////////////////////////
      fens[0] = e_fens[0];
      fens[1] = e_fens[1];
      fens[2] = e_fens[2];
      fens[3] = e_fens[4];
      _child[nc] = new GCELL_SOLID_T4 (fens);
      _child[nc]->_parent = this;
      this->gcell_group ()->add (_child[nc]);
      _child[nc]->_is_reflection = false;
      nc = 6; ///////////////////////////////////
      fens[0] = e_fens[5];
      fens[1] = e_fens[2];
      fens[2] = e_fens[1];
      fens[3] = e_fens[4];
      _child[nc] = new GCELL_SOLID_T4 (fens);
      _child[nc]->_parent = this;
      this->gcell_group ()->add (_child[nc]);
      _child[nc]->_is_reflection = true;
      nc = 7; ///////////////////////////////////
      fens[0] = e_fens[5];
      fens[1] = e_fens[4];
      fens[2] = e_fens[3];
      fens[3] = e_fens[2];
      _child[nc] = new GCELL_SOLID_T4 (fens);
      _child[nc]->_parent = this;
      this->gcell_group ()->add (_child[nc]);
      _child[nc]->_is_reflection = true;
    }
  }
}


void 
GCELL_SOLID_T4::detail_set (FEN *fen_of_bfun_to_refine,set <FEN*> *rf)
{
  // set <FEN *> result; result.clear();
  int local_index = _conn.local_index (fen_of_bfun_to_refine);
  if (local_index >= 0) {  // Insert nodes
    const sizet nfens = 4;
    // Insert nodes (they are just the nodes of the child adjacent
    // to the node whose bfun is to be refined
    for (sizet i = 0; i < nfens; i++) {
      if (i != (sizet) local_index) {
        // result.insert (_child[local_index]->conn()->fen(i));
        (*rf).insert (_child[local_index]->conn()->fen(i));
      }
    }
  }
  // return result;
}



void
GCELL_SOLID_T4::complete_refinement_set (FEN *fen_of_bfun_to_refine,set <FEN*> *rf )
{
  int local_index = _conn.local_index (fen_of_bfun_to_refine);
  if (local_index >= 0) {  // Insert nodes
    const sizet nfens = 4;
    // Insert nodes (they are just the nodes of the child adjacent
    // to the node whose bfun is to be refined
    for (sizet i = 0; i < nfens; i++) {
      (*rf).insert (_child[local_index]->conn()->fen(i));
      //gf->observe(_child[local_index]->conn()->fen(i));
    }
  }
}

void
GCELL_SOLID_T4::map_fen (FEN *fen, POINT *param_loc)
{
  if (fen->id() == _conn.fen(0)->id()) {
    *param_loc = 0;
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
    (*param_loc)(0) = 0; (*param_loc)(1) = 0; (*param_loc)(2) = +1;
    return;
  }
  throw EXCEPTION_BAD_VALUE();  // the node is unknown to this cell; raise hell
}

bool
GCELL_SOLID_T4::map_to_child (POINT &param_loc, GCELL **ret_child, POINT *child_param_loc)
{
  sizet j;
  GCELL_SOLID_T4 *child = 0;
  if (nchildren() == 0) return false;  // no children
  double xi = param_loc(0), eta = param_loc(1), theta = param_loc(2);
  if (xi > 0.5) {
    child = _child[1];
    j = 1;
  } else {
    if (eta > 0.5) {
      child = _child[2];
      j = 2;
    } else {
      if (theta > 0.5) {
        child = _child[3];
        j = 3;
      } else {
        // if greater than eq for plan (inside of plane) it is child 0
        if (0.5 - xi - eta - theta > 0) {
          child = _child[0];
          j = 0;
        } else {

        // if the program gets this far, then the child is one of the children
        // in the octahedron

          // check to see if point is inside (greater than 0) of plane 4,6,8,9
          if (0.5 - xi - (eta - xi) - theta > 0) {

            // check to see if this tetrahedra is a reflection of original, this will determine which
            // plane needs to be used for the determination of the child

            if (this->_is_reflection == false) {

              // check to see if point is inside of plane 4,5,7,9
              if (0.5 - (xi - eta) - eta - theta > 0) {
                // if inside of both planes 4,5,7,9 & 4,6,8,9 (and it is not children 0-3) then...
                child = _child[7];
                j = 7;
              } else {
                // if inside of planes 4,5,7,9 & outside 4,6,8,9 then...
                child = _child[6];
                j = 6;
              }

            } else { // reflection

              // the tetrahedra is a reflection so check to see if point is inside 5,6,7,8
              if (0.5 - xi - (eta - xi) - theta > 0) {
                // if inside 4,6,8,9 and inside 5,6,7,8
                child = _child[4];
                j = 4;
              } else {
                // if inside 4,6,8,9 and outside 5,6,7,8
                child = _child[7];
                j = 7;
              }

            } // reflection

          } else {  // outside of plane 4,6,8,9

            if (this->_is_reflection == false) {

              // check to see if point is inside of plane 4,5,7,9
              if (0.5 - (xi - eta) - eta - theta > 0) {
                // if outside of plane 4,5,7,9  and inside of 4,6,8,9 then...
                child = _child[4];
                j = 4;
              } else {
                // if outside of plane 4,5,7,9 and outside plane 4,6,8,9 then...
                child = _child[5];
                j = 5;
              }

            } else {  // reflection

              // the tetrahedra is a reflection so check to see if point is inside 5,6,7,8
              if (0.5 - xi - (eta - xi) - theta > 0) {
                // if outside 4,6,8,9 and inside 5,6,7,8
                child = _child[5];
                j = 5;
              } else {
                // if outside 4,6,8,9 and outside 5,6,7,8
                child = _child[6];
                j = 6;
              }

            } // reflection

          } // plane 4,6,8,9
        } // child 0
      } // child 3
    } // child 2
  } // child 1
  if ((j < 4) || (this->_is_reflection == false)) {
    (*child_param_loc)(0)
      = MPC[j].B[0][0] * xi + MPC[j].B[0][1] * eta + MPC[j].B[0][2] * theta + MPC[j].t[0];
    (*child_param_loc)(1)
      = MPC[j].B[1][0] * xi + MPC[j].B[1][1] * eta + MPC[j].B[1][2] * theta + MPC[j].t[1];
    (*child_param_loc)(2)
      = MPC[j].B[2][0] * xi + MPC[j].B[2][1] * eta + MPC[j].B[2][2] * theta + MPC[j].t[2];
  } else {
    j -= 4;
    (*child_param_loc)(0)
      = MPCR[j].B[0][0] * xi + MPCR[j].B[0][1] * eta + MPCR[j].B[0][2] * theta + MPCR[j].t[0];
    (*child_param_loc)(1)
      = MPCR[j].B[1][0] * xi + MPCR[j].B[1][1] * eta + MPCR[j].B[1][2] * theta + MPCR[j].t[1];
    (*child_param_loc)(2)
      = MPCR[j].B[2][0] * xi + MPCR[j].B[2][1] * eta + MPCR[j].B[2][2] * theta + MPCR[j].t[2];
  }
  *ret_child = child;
  return true;
}

bool
GCELL_SOLID_T4::map_to_parent (POINT &param_loc, GCELL **parent, POINT *parent_param_loc)
{
  if (_parent) {
    double xi = param_loc(0), eta = param_loc(1), theta = param_loc(2);
    for (sizet j = 0; j < GCELL_SOLID_T4::NCHILDREN; j++) {
      if (this == _parent->child(j)) {
        if ((j < 4) || (! _parent->_is_reflection)) {
          (*parent_param_loc)(0)
            = MCP[j].B[0][0] * xi + MCP[j].B[0][1] * eta + MCP[j].B[0][2] * theta + MCP[j].t[0];
          (*parent_param_loc)(1)
            = MCP[j].B[1][0] * xi + MCP[j].B[1][1] * eta + MCP[j].B[1][2] * theta + MCP[j].t[1];
          (*parent_param_loc)(2)
            = MCP[j].B[2][0] * xi + MCP[j].B[2][1] * eta + MCP[j].B[2][2] * theta + MCP[j].t[2];
        } else {
          j -= 4;
          (*parent_param_loc)(0)
            = MCPR[j].B[0][0] * xi + MCPR[j].B[0][1] * eta + MCPR[j].B[0][2] * theta + MCPR[j].t[0];
          (*parent_param_loc)(1)
            = MCPR[j].B[1][0] * xi + MCPR[j].B[1][1] * eta + MCPR[j].B[1][2] * theta + MCPR[j].t[1];
          (*parent_param_loc)(2)
            = MCPR[j].B[2][0] * xi + MCPR[j].B[2][1] * eta + MCPR[j].B[2][2] * theta + MCPR[j].t[2];
        }
       *parent = _parent;
        return true;
      }
    }
  }
  return false; // no parent
}

void
GCELL_SOLID_T4::debug_display ()
{
  cerr << "GCELL_SOLID_T4 " 
       << this->conn()->fen(0)->id() << " " << this->conn()->fen(0)->ref_loc() << " "
       << this->conn()->fen(1)->id() << " " << this->conn()->fen(1)->ref_loc() << " "
       << this->conn()->fen(2)->id() << " " << this->conn()->fen(2)->ref_loc() << " "
       << this->conn()->fen(3)->id() << " " << this->conn()->fen(3)->ref_loc() << " ";
  if (_parent) {
    cerr << endl;
    cerr << "parent= " 
       << this->_parent->conn()->fen(0)->id() << " "
       << this->_parent->conn()->fen(1)->id() << " "
       << this->_parent->conn()->fen(2)->id() << " "
       << this->_parent->conn()->fen(3)->id() << " ";
  }
    
  for (sizet j = 0; j < GCELL_SOLID_T4::NCHILDREN; j++) {
    if (_child[j]) {
      cerr << endl;
      cerr << "child[" << j << "]= " 
           << this->_child[j]->conn()->fen(0)->id() << " "
           << this->_child[j]->conn()->fen(1)->id() << " "
           << this->_child[j]->conn()->fen(2)->id() << " "
           << this->_child[j]->conn()->fen(3)->id() << " ";
    }
  }
}

#if 0
void
GCELL_SOLID_T4::check_is_inside_child (double xi, double eta, double theta, GCELL_SOLID_T4 *child,
                                       double child_xi, double child_eta, double child_theta)
{
  double N[4];
  N[0] = 1 - xi - eta - theta; N[1] = xi; N[2] = eta; N[3] = theta; 
  POINT loc;
  loc.assign (N[0], _conn.fen(0)->ref_loc(),
              N[1], _conn.fen(1)->ref_loc(),
              N[2], _conn.fen(2)->ref_loc(),
              N[3], _conn.fen(3)->ref_loc());
  N[0] = 1 - child_xi - child_eta - child_theta; N[1] = child_xi; N[2] = child_eta; N[3] = child_theta; 
  G3D_loc_t pt;
  pt.x = loc(0);
  pt.y = loc(1);
  pt.z = loc(2);
  POINT child_loc;
  child_loc.assign (N[0], child->_conn.fen(0)->ref_loc(),
                    N[1], child->_conn.fen(1)->ref_loc(),
                    N[2], child->_conn.fen(2)->ref_loc(),
                    N[3], child->_conn.fen(3)->ref_loc());
  G3D_tet_t t;
  t.loc[0].x = child->_conn.fen(0)->ref_loc()(0);
  t.loc[0].y = child->_conn.fen(0)->ref_loc()(1);
  t.loc[0].z = child->_conn.fen(0)->ref_loc()(2);
  t.loc[1].x = child->_conn.fen(1)->ref_loc()(0);
  t.loc[1].y = child->_conn.fen(1)->ref_loc()(1);
  t.loc[1].z = child->_conn.fen(1)->ref_loc()(2);
  t.loc[2].x = child->_conn.fen(2)->ref_loc()(0);
  t.loc[2].y = child->_conn.fen(2)->ref_loc()(1);
  t.loc[2].z = child->_conn.fen(2)->ref_loc()(2);
  t.loc[3].x = child->_conn.fen(3)->ref_loc()(0);
  t.loc[3].y = child->_conn.fen(3)->ref_loc()(1);
  t.loc[3].z = child->_conn.fen(3)->ref_loc()(2);
  if (!G3D_point_is_inside_tet (&pt, &t, 1e-6)) {
    cerr << "Error mapping ***** in" << __FUNCTION__ << "line " << __LINE__ << endl;
    cerr << "Error mapping ***** " << endl;
    cerr << "  parent " << loc << endl;
    cerr << "  child " << child_loc << endl;
  } else {
    cerr << "Passed in " << __FUNCTION__ << endl;
  }
  POINT diff;
  diff = child_loc - loc;
  if (diff.l2_norm() > 1e-6) {
    cerr << "Error mapping ***** in " << __FUNCTION__ << "line " << __LINE__ << endl;
    cerr << "  parent " << loc << "[" << xi << " " << eta << " " << theta << "]" << endl;
    cerr << "  child " << child_loc << "[" << child_xi << " " << child_eta << " " << child_theta << "]" << endl;
    for (sizet j = 0; j < 8; j++) {
      if (child == _child[j]) cerr << "child #" << j << endl;
    }
    if (this->_is_reflection) cerr << "this is a reflection" << endl;
    if (child->_is_reflection) cerr << "child is a reflection" << endl;
    this->debug_display();
    child->debug_display();
  } else {
    cerr << "Passed in " << __FUNCTION__ << endl;
  }
}

void
GCELL_SOLID_T4::check_is_inside_parent (double xi, double eta, double theta, GCELL_SOLID_T4 *parent,
                                       double parent_xi, double parent_eta, double parent_theta)
{
  double N[4];
  N[0] = 1 - xi - eta - theta; N[1] = xi; N[2] = eta; N[3] = theta; 
  POINT loc;
  loc.assign (N[0], _conn.fen(0)->ref_loc(),
              N[1], _conn.fen(1)->ref_loc(),
              N[2], _conn.fen(2)->ref_loc(),
              N[3], _conn.fen(3)->ref_loc());
  N[0] = 1 - parent_xi - parent_eta - parent_theta; N[1] = parent_xi; N[2] = parent_eta; N[3] = parent_theta; 
  POINT parent_loc;
  parent_loc.assign (N[0], parent->_conn.fen(0)->ref_loc(),
                    N[1], parent->_conn.fen(1)->ref_loc(),
                    N[2], parent->_conn.fen(2)->ref_loc(),
                    N[3], parent->_conn.fen(3)->ref_loc());
  POINT diff = loc - parent_loc;
  if (diff.l2_norm() > 1e-6) {
    cerr << "Error mapping ***** in" << __FUNCTION__ << endl;
    cerr << "  parent " << parent_loc << endl;
    cerr << "  child " << loc << endl;
  } else {
    cerr << "Passed in " << __FUNCTION__ << endl;
  }
}
#endif
