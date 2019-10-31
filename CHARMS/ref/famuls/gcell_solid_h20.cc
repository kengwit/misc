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
#include "gcell_solid_h20.h"
#include "field.h"
#include "evalpt.h"
#include "ref_ctx.h"

const char GCELL_SOLID_H20::TYPE_NAME[] = "solid_h20";

const double GCELL_SOLID_H20::_child_map[NCHILDREN][2][3] = {
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
#include "g3dhex20.h"
}

GCELL_SOLID_H20::GCELL_SOLID_H20 (vector <FEN *> fens) : GCELL ()
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
  return (new GCELL_SOLID_H20 (fens));
}

bool
GCELL_SOLID_H20::register_make_func ()
{
  return GCELL::register_make_func (make, string(TYPE_NAME), string ("default"));
}

void
GCELL_SOLID_H20::eval_ns (GCELL_SOLID_H20 *child, EVALPT *evalpt, 
                         double xi, double eta, double theta, double dxiparent_dxichild)
{
  const sizet nbfuns = evalpt->nbfuns();
  if (child != 0) {
    for (sizet j = 0; j < GCELL_SOLID_H20::NCHILDREN; j++) {
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
  const sizet nfens = 20;
  double N[nfens];
  G3D_hex20_N (N, xi, eta, theta);
  double N_xi[nfens], N_eta[nfens], N_theta[nfens];
  G3D_hex20_N_xi    (   N_xi, xi, eta, theta);
  G3D_hex20_N_eta   (  N_eta, xi, eta, theta);
  G3D_hex20_N_theta (N_theta, xi, eta, theta);
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
GCELL_SOLID_H20::eval_bfun_set (EVALPT *evalpt)
{
  POINT param_loc = evalpt->param_loc();
  eval_ns (0, evalpt, param_loc(0), param_loc(1), param_loc(2), 1.0);
  return true;
}

void
GCELL_SOLID_H20::divide (REF_CTX *ref_ctx)
{
  CHECK(0, EXCEPTION_NOT_IMPLEMENTED,;);
}


void
GCELL_SOLID_H20::detail_set (FEN *fen_of_bfun_to_refine,set <FEN *> *rf)
{
  CHECK(0, EXCEPTION_NOT_IMPLEMENTED,;);
}



void
GCELL_SOLID_H20::complete_refinement_set (FEN *fen_of_bfun_to_refine,set <FEN *> *rf)
{
  CHECK(0, EXCEPTION_NOT_IMPLEMENTED,;);
}


extern "C" {
#include "g3dhex20.h"
}
void
GCELL_SOLID_H20::map_fen (FEN *fen, POINT *param_loc)
{
  for (sizet j = 0; j < _conn.nfens (); j++) {
    if (fen->id() == _conn.fen(j)->id()) {
      (*param_loc)(0) = g3dhex20_vertex_param_coord[j][0];
      (*param_loc)(1) = g3dhex20_vertex_param_coord[j][1];
      (*param_loc)(2) = g3dhex20_vertex_param_coord[j][2];
      return;
    }
  }
  throw EXCEPTION_BAD_VALUE();  // the node is unknown to this cell; raise hell
}

bool
GCELL_SOLID_H20::map_to_child (POINT &param_loc, GCELL **child, POINT *child_param_loc)
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
GCELL_SOLID_H20::map_to_parent (POINT &param_loc, GCELL **parent, POINT *parent_param_loc)
{
  CHECK(0, EXCEPTION_NOT_IMPLEMENTED,;);
  return false; // no parent
}

void
GCELL_SOLID_H20::debug_display ()
{
  cerr << "GCELL_SOLID_H20 "
       << this->conn()->fen(0)->id() << " "
       << this->conn()->fen(1)->id() << " "
       << this->conn()->fen(2)->id() << " "
       << this->conn()->fen(3)->id() << " "
       << this->conn()->fen(4)->id() << " "
       << this->conn()->fen(5)->id() << " "
       << this->conn()->fen(6)->id() << " "
       << this->conn()->fen(7)->id() << " "
       << this->conn()->fen(8)->id() << " "
       << this->conn()->fen(9)->id() << " "
       << this->conn()->fen(10)->id() << " "
       << this->conn()->fen(11)->id() << " "
       << this->conn()->fen(12)->id() << " "
       << this->conn()->fen(13)->id() << " "
       << this->conn()->fen(14)->id() << " "
       << this->conn()->fen(15)->id() << " "
       << this->conn()->fen(16)->id() << " "
       << this->conn()->fen(17)->id() << " "
       << this->conn()->fen(18)->id() << " "
       << this->conn()->fen(19)->id() << " "
       << this->conn()->fen(20)->id();
}

