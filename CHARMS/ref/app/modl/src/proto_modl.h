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
#ifndef PROTO_MODL_H
# define PROTO_MODL_H

#include <list>
#include "algo.h"
#include "gmesh.h"
#include "ecell_modl.h"
#include "evs_ooofs.h"
#include "mat_mgr.h"
#include "mat_map.h"
#include "mat_elasticity.h"

/**
  Prototypical modal equation: K u = omega^2 M u

 */
class PROTO_MODL {

 public: // object functions ////////////////////////////////////////

  PROTO_MODL (ALGO *algo,
                          DB *db,
                          EVS_OOOFS *evs, 
                          FIELD_VECTOR *geometry, 
                          FIELD_VECTOR *u 
                          );
  ~PROTO_MODL () ;
  /**
     Assemble the stiffness matrix K (see the class comment).
  */
  bool assemble_stiffness_matrix (double pressure, double shift);
  bool assemble_stiffness_matrix ( double shift);
  bool assemble_mass_matrix(); 
  double total_mass(); 
  double total_strain_energy_in_gcell_group(FIELD_VECTOR *u, GCELL_GROUP *gg);
  /**
     Assemble the load (rhs) terms corresponding to non-zero prescribed
     boundary displacements, applied body loads, and applied tractions.
  */
  bool assemble_load_terms ();
  /**
     Solve for the unknown field u.  (Both the stiffness matrix
     and the load terms must have been assembled prior to this call.)
  */
  bool solve (double shift);
  /**
     Get the EVS (Eigen Value Solver) that is being
     used by the protocol to solve the diffusion equation.
  */
  EVS_OOOFS *evs () { return _evs; }
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
  bool assemble_mass ();
  
  
 private: // object data /////////////////////////////////////////////

  ALGO                           *_algo;
  DB                             *_db;
  EVS_OOOFS                      *_evs;
  FIELD_VECTOR                   *_geometry;
  list <ECELL_MODL *>             _ecells; // responsible for allocation/deallocation
  FIELD_VECTOR                   *_u;
  MAT_MAP <MAT_ELASTICITY >      *_mat_map;

};

#endif
