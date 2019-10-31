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
#include "gcell_solid_t10.h"
#include "field.h"
#include "evalpt.h"
#include "ref_ctx.h"

const char GCELL_SOLID_T10::TYPE_NAME[] = "solid_t10";

const GCELL_SOLID_T10::map_struct GCELL_SOLID_T10::map_child_to_parent[NCHILDREN] = {
{{{0.5,0,0},{0,0.5,0},{0,0,0.5}},{0,0,0}}, // child 0
{{{0.5,0,0},{0,0.5,0},{0,0,0.5}},{0.5,0,0}}, // child 1
{{{0.5,0,0},{0,0.5,0},{0,0,0.5}},{0,0.5,0}}, // child 2
{{{0.5,0,0},{0,0.5,0},{0,0,0.5}},{0,0,0.5}}, // child 3
{{{-0.5,-0.5,0},{0,0.5,0},{0,0,-0.5}},{0.5,0,0.5}}, // child 4
{{{0,0,-0.5},{0.5,0,0.5},{-0.5,-0.5,0}},{0.5,0,0.5}}, // child 5
{{{0.5,0.5,0},{-0.5,0,0},{0,0,0.5}},{0,0.5,0}}, // child 6
{{{0,0.5,0},{0,-0.5,-0.5},{-0.5,-0.5,0}},{0,0.5,0.5}}, // child 7
};

const GCELL_SOLID_T10::map_struct GCELL_SOLID_T10::map_child_to_parent_reflected[4] = {
{{{-0.5,0,-0.5},{0,0,0.5},{0.5,0.5,0}},{0.5,0,0}}, // child 4
{{{0,-0.5,0},{0.5,0.5,0},{0,0,0.5}},{0.5,0,0}}, // child 5
{{{0,0.5,0.5},{0,0,-0.5},{-0.5,-0.5,0}},{0,0.5,0.5}}, // child 6
{{{0.5,0,0},{-0.5,-0.5,0},{0,0,-0.5}},{0,0.5,0.5}}, // child 7
};

const GCELL_SOLID_T10::map_struct GCELL_SOLID_T10::map_parent_to_child[NCHILDREN] = {
{{{2,-0,-0},{0,2,-0},{0,0,2}},{-0,-0,-0}}, // child 0
{{{2,-0,-0},{0,2,-0},{0,0,2}},{-1,-0,-0}}, // child 1
{{{2,-0,-0},{0,2,-0},{0,0,2}},{-0,-1,-0}}, // child 2
{{{2,-0,-0},{0,2,-0},{0,0,2}},{-0,-0,-1}}, // child 3
{{{-2,-2,-0},{0,2,0},{0,0,-2}},{1,-0,1}}, // child 4
{{{2,2,0},{-2,-2,-2},{-2,0,0}},{-1,2,1}}, // child 5
{{{0,-2,-0},{2,2,-0},{0,0,2}},{1,-1,-0}}, // child 6
{{{-2,0,-2},{2,0,0},{-2,-2,0}},{1,-0,1}}, // child 7
};

const GCELL_SOLID_T10::map_struct GCELL_SOLID_T10::map_parent_to_child_reflected[4] = {
{{{-2,-2,0},{2,2,2},{0,2,0}},{1,-1,-0}}, // child 4
{{{2,2,-0},{-2,0,0},{0,0,2}},{-1,1,-0}}, // child 5
{{{-2,-2,-2},{2,2,0},{0,-2,0}},{2,-1,1}}, // child 6
{{{2,0,0},{-2,-2,-0},{0,0,-2}},{-0,1,1}}, // child 7
};

const double GCELL_SOLID_T10::detfen_param_c[25][3] = {
  {0.25,0,0}, // 10
  {0.75,0,0}, // 11
  {0.75,0.25,0}, // 12
  {0.25,0.75,0}, // 13
  {0,0.75,0}, // 14
  {0,0.25,0}, // 15
  {0,0,0.25}, // 16
  {0,0,0.75}, // 17
  {0.75,0,0.25}, // 18
  {0.25,0,0.75}, // 19
  {0,0.75,0.25}, // 20
  {0,0.25,0.75}, // 21
  {0.25,0.25,0}, // 22
  {0.25,0.50,0}, // 23
  {0.50,0.25,0}, // 24
  {0.25,0,0.25}, // 25
  {0.50,0,0.25}, // 26
  {0.25,0,0.50}, // 27
  {0.50,0.25,0.25}, // 28
  {0.25,0.50,0.25}, // 29
  {0.25,0.25,0.50}, // 30
  {0,0.50,0.25}, // 31
  {0,0.25,0.25}, // 32
  {0,0.25,0.50}, // 33
  {0.25,0.25,0.25} // 34
  };

#define MCP map_child_to_parent
#define MCPR map_child_to_parent_reflected
#define MPC map_parent_to_child
#define MPCR map_parent_to_child_reflected

// The first 10 nodes are those that result from the addition
// of six refinement nodes to the four of the t4 gcell.
// Then, on the edges we have nodes 10 thru 21
// 0---10---4---11---1
// 1---12---5---13---2
// 2---14---6---15---0
// 0---16---7---17---3
// 1---18---8---19---3
// 2---20---9---21---3
// and on the faces
// 
// and on the body diagonal: 34
const int GCELL_SOLID_T10::_child_map[NCHILDREN][10] = {
  {0,4,6,7,10,22,15,16,25,32}, // child 0
  {4,1,5,8,11,12,24,26,18,28}, // child 1
  {6,5,2,9,23,13,14,31,29,20}, // child 2
  {7,8,9,3,27,30,33,17,19,21}, // child 3
  {8,7,9,4,27,33,30,26,25,34}, // child 4
  {8,5,4,9,28,24,26,30,29,34}, // child 5
  {6,4,5,9,22,24,23,31,34,29}, // child 6
  {9,6,4,7,31,22,34,33,32,25}  // child 7
};

// 
const int GCELL_SOLID_T10::_child_map_reflected[NCHILDREN][10] = {
  {0,4,6,7,10,22,15,16,25,32}, // child 0, same as above
  {4,1,5,8,11,12,24,26,18,28}, // child 1, same as above
  {6,5,2,9,23,13,14,31,29,20}, // child 2, same as above
  {7,8,9,3,27,30,33,17,19,21}, // child 3, same as above
  {4,7,8,6,25,27,26,22,32,34}, // child 4
  {4,5,6,8,24,23,22,26,28,34}, // child 5
  {9,6,5,8,31,23,29,30,34,28}, // child 6
  {9,8,7,6,30,27,33,31,34,32}  // child 7
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


GCELL_SOLID_T10::GCELL_SOLID_T10 (vector <FEN *> fens) : GCELL ()
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
  return (new GCELL_SOLID_T10 (fens));
}

bool
GCELL_SOLID_T10::register_make_func ()
{
  return GCELL::register_make_func (make, string(TYPE_NAME), string ("default"));
}

void
GCELL_SOLID_T10::eval_ns (GCELL_SOLID_T10 *child, EVALPT *evalpt,
                         double xi, double eta, double theta, double B[3][3])
{
  const sizet nbfuns = evalpt->nbfuns ();
  // RAW: Should be told what needs to be evaluated (connectivity, basis function value, basis
  //      function value + derivatives) in order to be able to optimize the effort
  if (child != 0) {
    for (sizet j = 0; j < GCELL_SOLID_T10::NCHILDREN; j++) {
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
  const sizet nfens = 10;
  double N[nfens]; // BELOW (the basis functions)
  double upsilon = 1 - xi - eta - theta;
  G3D_tet10_N (N, xi, eta, theta, upsilon);
  double N_xi[nfens], N_eta[nfens], N_theta[nfens];
  G3D_tet10_N_r (N_xi, xi, eta, theta, upsilon);
  G3D_tet10_N_s (N_eta, xi, eta, theta, upsilon);
  G3D_tet10_N_t (N_theta, xi, eta, theta, upsilon);
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
GCELL_SOLID_T10::eval_bfun_set (EVALPT *evalpt)
{
  double B[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
  POINT param_loc = evalpt->param_loc();
  eval_ns (0, evalpt, param_loc(0), param_loc(1), param_loc(2), B);

  return true;
}

void
GCELL_SOLID_T10::divide (REF_CTX *ref_ctx)
{
  if (this->nchildren () == 0) { // Not refined yet
    const sizet nfens = 10;
    FEN *r_fens[35];
    // Get the refinement nodes 
    int rc = 0;
    // The first ten are at the original nodes of the parent
    {
      CONN_POINT_1 vertex_conn;
      vector <POINT> ref_locs(1);
      for (sizet i = 0; i < nfens; i++) {
        vertex_conn.fen(0) = _conn.fen(i);
        ref_locs[0] = _conn.fen(i)->ref_loc();
        r_fens[rc] = ref_ctx->get_ref_fen (&vertex_conn, 0, ref_locs); rc++;
      }
    }
    // Get the refinement nodes for the edges
    {
      CONN_LINE_3 edge_conn;
      vector <POINT> ref_locs(2);
      const double refnrs[2] = {-0.5,+0.5};
      for (int i = 0; i < 6; i++) {
        edge_conn.fen(0) = _conn.fen(g3dtet10_vertex_on_edge[i][0]);
        edge_conn.fen(1) = _conn.fen(g3dtet10_vertex_on_edge[i][2]);
        edge_conn.fen(2) = _conn.fen(g3dtet10_vertex_on_edge[i][1]); // mid-point
        for (int k = 0; k < 2; k++) {
          double n[3];
          G3D_line3_N (n, refnrs[k]);
          ref_locs[k].assign (n[0], edge_conn.fen(0)->ref_loc(),
                              n[1], edge_conn.fen(2)->ref_loc(),
                              n[2], edge_conn.fen(1)->ref_loc());
        }
        for (int k = 0; k < 2; k++) {
          r_fens[rc] = ref_ctx->get_ref_fen (&edge_conn, k, ref_locs); rc++;
        }
      }
    }
    // Get the refinement nodes for the faces
    {
      CONN_SURF_6 face_conn;
      vector <POINT> ref_locs(3);
      const double refnrs[3][2] = { // parametric coordinates of face refinement nodes
        {0.25, 0.25},
        {0.50, 0.25},
        {0.25, 0.50}
      };
      for (int i = 0; i < 4; i++) {
        face_conn.fen(0) = _conn.fen(g3dtet10_vertex_on_face[i][0]);
        face_conn.fen(1) = _conn.fen(g3dtet10_vertex_on_face[i][2]);
        face_conn.fen(2) = _conn.fen(g3dtet10_vertex_on_face[i][4]); 
        face_conn.fen(3) = _conn.fen(g3dtet10_vertex_on_face[i][1]); 
        face_conn.fen(4) = _conn.fen(g3dtet10_vertex_on_face[i][3]); 
        face_conn.fen(5) = _conn.fen(g3dtet10_vertex_on_face[i][5]); 
        for (int k = 0; k < 3; k++) {
          double n[6];
          G3D_tri6_N (n, refnrs[k][0], refnrs[k][1]);
          ref_locs[k] = 0;
          for (int j = 0; j < 6; j++)
            ref_locs[k].add (n[j], face_conn.fen(j)->ref_loc());
        }
        for (int k = 0; k < 3; k++) {
          r_fens[rc] = ref_ctx->get_ref_fen (&face_conn, k, ref_locs); rc++;
        }
      }
    }
    // Get the refinement node for the interior
    {
      vector <POINT> ref_locs(1);
      double n[10];
      G3D_tet10_N (n, 0.25, 0.25, 0.25, 0.25);
      ref_locs[0] = 0;
      for (int j = 0; j < 10; j++) {
        ref_locs[0].add (n[j], _conn.fen(j)->ref_loc());
      }
      r_fens[rc] = ref_ctx->get_ref_fen (&_conn, 0, ref_locs); rc++;
    }
    // Now generate the new gcells
    vector <FEN *> fens(10);
    // First the four corner children
    for (int nc = 0; nc < 4; nc++) {
      for (int j = 0; j < 10; j++) { fens[j] = r_fens[_child_map[nc][j]]; }
      _child[nc] = new GCELL_SOLID_T10 (fens);
      _child[nc]->_parent = this;
      this->gcell_group ()->add (_child[nc]);
      _child[nc]->_is_reflection = this->_is_reflection;
      //_child[nc]->debug_display();
    }
    // The mapping of the four octahedral children depends on whether
    // the parent is a reflection or not.
    const int is_reflection[4] = {false, false, true, true};
    if (this->_is_reflection) {
      for (int nc = 4; nc < 8; nc++) {
        for (int j = 0; j < 10; j++) { fens[j] = r_fens[_child_map_reflected[nc][j]]; }
        _child[nc] = new GCELL_SOLID_T10 (fens);
        _child[nc]->_parent = this;
        this->gcell_group ()->add (_child[nc]);
        _child[nc]->_is_reflection = is_reflection[nc-4];
        //_child[nc]->debug_display();
      }
    } else {
      for (int nc = 4; nc < 8; nc++) {
        for (int j = 0; j < 10; j++) { fens[j] = r_fens[_child_map[nc][j]]; }
        _child[nc] = new GCELL_SOLID_T10 (fens);
        _child[nc]->_parent = this;
        this->gcell_group ()->add (_child[nc]);
        _child[nc]->_is_reflection = is_reflection[nc-4];
        //_child[nc]->debug_display();
      }
    }
  }
}


void
GCELL_SOLID_T10::detail_set (FEN *fen_of_bfun_to_refine, set <FEN *> *rf)
{
  int local_index = _conn.local_index (fen_of_bfun_to_refine);
  if (local_index >= 0) {  // Insert nodes
#undef A
# define A(whichfen) (*rf).insert (r_fens[whichfen])
#undef A13
# define A13(i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12) { \
   A(i0);A(i1);A(i2);A(i3);A(i4);A(i5);A(i6);A(i7);A(i8);A(i9);A(i10);A(i11);A(i12); \
  }
#undef A9
# define A9(i0,i1,i2,i3,i4,i5,i6,i7,i8) { \
   A(i0);A(i1);A(i2);A(i3);A(i4);A(i5);A(i6);A(i7);A(i8); \
  }
    // Collect refinement nodes
    FEN *r_fens[35];
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 10; k++) {
        int m = _child_map[j][k];
        r_fens[m] = this->_child[j]->_conn.fen(k);
      }
    }
    r_fens[34] = this->_child[4]->_conn.fen(9); // always the last node in child 4:
                                                //see _child_map and _child_map_reflected
    switch (local_index) {
    case 0: // corner
      A13(10,15,16,17,27,26,11,24,23,14,31,33,34);
      break;
    case 1:
      A13(11,12,18,10,22,23,13,29,30,19,27,25,34);
      break;
    case 2:
      A13(13,14,20,12,24,22,15,32,33,21,30,28,34);
      break;
    case 3:
      A13(17,19,21,16,25,26,18,28,29,20,31,32,34);
      break;
    case 4:
      A9(11,10,22,23,24,25,26,27,34);
      break;
    case 5:
      A9(12,13,22,23,24,28,29,30,34);
      break;
    case 6:
      A9(14,15,22,23,24,31,32,33,34);
      break;
    case 7:
      A9(16,17,25,26,27,31,32,33,34);
      break;
    case 8:
      A9(18,19,25,26,27,28,29,30,34);
      break;
    case 9:
      A9(20,21,28,29,30,31,32,33,34);
      break;
    }
  }
}


void
GCELL_SOLID_T10::complete_refinement_set (FEN *fen_of_bfun_to_refine, set <FEN *> *rf)
{
  int local_index = _conn.local_index (fen_of_bfun_to_refine);
  if (local_index >= 0) {  // Insert nodes
    detail_set (fen_of_bfun_to_refine, rf);
    FEN *r_fens[35];
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 10; k++) {
        int m = _child_map[j][k];
        r_fens[m] = this->_child[j]->_conn.fen(k);
      }
    }
    r_fens[34] = this->_child[4]->_conn.fen(9); // always the last node in child 4:
                                                //see _child_map and _child_map_reflected
    (*rf).insert(r_fens[local_index]);
  }
}

void
GCELL_SOLID_T10::map_fen (FEN *fen, POINT *param_loc)
{
  for (int j = 0; j < 10; j++) {
    if (fen->id() == _conn.fen(j)->id()) {
      (*param_loc)(0) = g3dtet10_vertex_param_coords[j][0];
      (*param_loc)(1) = g3dtet10_vertex_param_coords[j][1];
      (*param_loc)(2) = g3dtet10_vertex_param_coords[j][2];
      return;
    }
  }
  throw EXCEPTION_BAD_VALUE();  // the node is unknown to this cell; raise hell
}

bool
GCELL_SOLID_T10::map_to_child (POINT &param_loc, GCELL **child, POINT *child_param_loc)
{
  sizet j;
  if (nchildren() == 0) return false;  // no children
  double xi = param_loc(0), eta = param_loc(1), theta = param_loc(2);
  if (xi > 0.5) {
    *child = _child[1];
    j = 1;
  } else {
    if (eta > 0.5) {
      *child = _child[2];
      j = 2;
    } else {
      if (theta > 0.5) {
        *child = _child[3];
        j = 3;
      } else {
        // if greater than eq for plan (inside of plane) it is child 0
        if (0.5 - xi - eta - theta > 0) {
          *child = _child[0];
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
                *child = _child[7];
                j = 7;
              } else {
                // if inside of planes 4,5,7,9 & outside 4,6,8,9 then...
                *child = _child[6];
                j = 6;
              }

            } else { // reflection

              // the tetrahedra is a reflection so check to see if point is inside 5,6,7,8
              if (0.5 - xi - (eta - xi) - theta > 0) {
                // if inside 4,6,8,9 and inside 5,6,7,8
                *child = _child[4];
                j = 4;
              } else {
                // if inside 4,6,8,9 and outside 5,6,7,8
                *child = _child[7];
                j = 7;
              }

            } // reflection

          } else {  // outside of plane 4,6,8,9

            if (this->_is_reflection == false) {

              // check to see if point is inside of plane 4,5,7,9
              if (0.5 - (xi - eta) - eta - theta > 0) {
                // if outside of plane 4,5,7,9  and inside of 4,6,8,9 then...
                *child = _child[4];
                j = 4;
              } else {
                // if outside of plane 4,5,7,9 and outside plane 4,6,8,9 then...
                *child = _child[5];
                j = 5;
              }

            } else {  // reflection

              // the tetrahedra is a reflection so check to see if point is inside 5,6,7,8
              if (0.5 - xi - (eta - xi) - theta > 0) {
                // if outside 4,6,8,9 and inside 5,6,7,8
                *child = _child[5];
                j = 5;
              } else {
                // if outside 4,6,8,9 and outside 5,6,7,8
                *child = _child[6];
                j = 6;
              }

            } // reflection

          } // plane 4,6,8,9
        } // child 0
      } // child 3
    } // child 2
  } // child 1
  if ((j < 4) || (this->_is_reflection == false)) {
    (*child_param_loc)(0) = MPC[j].B[0][0] * xi + MPC[j].B[0][1] * eta + MPC[j].B[0][2] * theta + MPC[j].t[0];
    (*child_param_loc)(1) = MPC[j].B[1][0] * xi + MPC[j].B[1][1] * eta + MPC[j].B[1][2] * theta + MPC[j].t[1];
    (*child_param_loc)(2) = MPC[j].B[2][0] * xi + MPC[j].B[2][1] * eta + MPC[j].B[2][2] * theta + MPC[j].t[2];
  } else {
    j -= 4;
    (*child_param_loc)(0) = MPCR[j].B[0][0] * xi + MPCR[j].B[0][1] * eta + MPCR[j].B[0][2] * theta + MPCR[j].t[0];
    (*child_param_loc)(1) = MPCR[j].B[1][0] * xi + MPCR[j].B[1][1] * eta + MPCR[j].B[1][2] * theta + MPCR[j].t[1];
    (*child_param_loc)(2) = MPCR[j].B[2][0] * xi + MPCR[j].B[2][1] * eta + MPCR[j].B[2][2] * theta + MPCR[j].t[2];
  }
  return true;
}

bool
GCELL_SOLID_T10::map_to_parent (POINT &param_loc, GCELL **parent, POINT *parent_param_loc)
{
  if (_parent) {
    double xi = param_loc(0), eta = param_loc(1), theta = param_loc(2);
    for (sizet j = 0; j < GCELL_SOLID_T10::NCHILDREN; j++) {
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
GCELL_SOLID_T10::debug_display ()
{
  cerr << "GCELL_SOLID_T10 " 
       << this->conn()->fen(0)->id() << " "
       << this->conn()->fen(1)->id() << " "
       << this->conn()->fen(2)->id() << " "
       << this->conn()->fen(3)->id() << " "
       << this->conn()->fen(4)->id() << " "
       << this->conn()->fen(5)->id() << " "
       << this->conn()->fen(6)->id() << " "
       << this->conn()->fen(7)->id() << " "
       << this->conn()->fen(8)->id() << " "
       << this->conn()->fen(9)->id() << " "; 
      cerr << endl;
  if (_parent) {
    cerr << "parent= " 
       << this->_parent->conn()->fen(0)->id() << " "
       << this->_parent->conn()->fen(1)->id() << " "
       << this->_parent->conn()->fen(2)->id() << " "
       << this->_parent->conn()->fen(3)->id() << " "
       << this->_parent->conn()->fen(4)->id() << " "
       << this->_parent->conn()->fen(5)->id() << " "
       << this->_parent->conn()->fen(6)->id() << " "
       << this->_parent->conn()->fen(7)->id() << " "
       << this->_parent->conn()->fen(8)->id() << " "
       << this->_parent->conn()->fen(9)->id() << " ";
    cerr << endl;
  }
    
  for (sizet j = 0; j < GCELL_SOLID_T10::NCHILDREN; j++) {
    if (_child[j]) {
      cerr << "child[" << j << "]= " 
           << this->_child[j]->conn()->fen(0)->id() << " "
           << this->_child[j]->conn()->fen(1)->id() << " "
           << this->_child[j]->conn()->fen(2)->id() << " "
           << this->_child[j]->conn()->fen(3)->id() << " "
           << this->_child[j]->conn()->fen(4)->id() << " "
           << this->_child[j]->conn()->fen(5)->id() << " "
           << this->_child[j]->conn()->fen(6)->id() << " "
           << this->_child[j]->conn()->fen(7)->id() << " "
           << this->_child[j]->conn()->fen(8)->id() << " "
           << this->_child[j]->conn()->fen(9)->id() << " ";
      cerr << endl;
    }
  }
}

