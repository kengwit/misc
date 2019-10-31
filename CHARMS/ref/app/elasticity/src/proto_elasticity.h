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
#ifndef PROTO_ELASTICITY_H
# define PROTO_ELASTICITY_H

#include <list>
#include "algo.h"
#include "gmesh.h"
#include "ecell_elasticity.h"
#include "select_les.h"
using ELASTICITY::LES;
#include "mat_mgr.h"
#include "mat_map.h"
#include "mat_elasticity.h"
#include "load_on_gcell.h"

/**
  Prototypical elasticity equation: K u = F

  K=stiffness matrix
  u=displacement vector
  F=load vector

 */
class PROTO_ELASTICITY {

 public: // object functions ////////////////////////////////////////

  PROTO_ELASTICITY (ALGO *algo,
                          DB *db,
                          LES *les, 
                          FIELD_VECTOR *geometry, 
                          FIELD_VECTOR *u 
                          );
  ~PROTO_ELASTICITY () ;
  /**
     Assemble the stiffness matrix K (see the class comment).
  */
  bool assemble_stiffness_matrix ();
  /**
     Assemble the load (rhs) terms corresponding to non-zero prescribed
     boundary displacements, applied body loads, and applied tractions.
  */
  bool assemble_load_terms ();
  /**
     Solve for the unknown field u.  (Both the stiffness matrix
     and the load terms must have been assembled prior to this call.)
  */
  bool solve ();
  /**
     Get the LES (Linear Equation Solver) that is being
     used by the protocol to solve the diffusion equation.
  */
  LES *les () { return _les; }
  /**
     Get the primary field that this protocol is used
     to solve for.
  */
  FIELD_VECTOR *u () { return _u; }
  /**
     Get the geometry field that this protocol is supposed
     to use in the solution.
  */
  FIELD_VECTOR *geometry () { return _geometry; }
  /**
     Get the associated database.
  */
  DB *db () { return _db; }
  /**
     Material.
  */
  MAT_ELASTICITY *mat (string mattype, string matname) ;
  /**
     Load.
  */
  LOAD_ON_GCELL *load (string loadtype, string loadname) ;
  
 private: // object data /////////////////////////////////////////////

  ALGO                           *_algo;
  DB                             *_db;
  LES                            *_les;
  FIELD_VECTOR                   *_geometry;
  list <ECELL_ELASTICITY *>       _ecells; // responsible for allocation/deallocation
  FIELD_VECTOR                   *_u;
  MAT_MAP <MAT_ELASTICITY >      *_mat_map;

};

#endif
