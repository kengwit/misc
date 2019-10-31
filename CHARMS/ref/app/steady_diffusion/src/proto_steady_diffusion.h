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
#ifndef PROTO_STEADY_DIFFUSION_H
# define PROTO_STEADY_DIFFUSION_H

#include <list>
#include "algo.h"
#include "gmesh.h"
#include "ecell_steady_diffusion.h"
#include "select_les.h"
#include "mat_mgr.h"
#include "mat_map.h"
#include "mat_diffusion_iso.h"

/**
  Prototypical heat equation: K(phi) phi = F(phi)

  d/dx_i (k_ij dphi/dx_j) + Q = 0

  C = specific heat, k_ij = cartesian components of the conductivity
  tensor (symmetric), Q = internal heat generation
  
  with bc's

  phi = f_phi (x) on Gamma_phi

  -(k_ij dphi/dx_j) n_i = q_a + q_c + q_r = f_q (x) on Gamma_q

  q_a = applied flux
  q_c = convective flux, q_c = h_c(x, phi) (phi - phi_c)
  q_r = radiative flux, q_r = h_r(x, phi) (phi - phi_r)

 */
class PROTO_STEADY_DIFFUSION {

 public: // object functions ////////////////////////////////////////

  PROTO_STEADY_DIFFUSION (ALGO *algo,
                          STEADY_DIFFUSION::LES *les, 
                          FIELD_VECTOR *geometry, 
                          FIELD_SCALAR *phi 
                          );
  /**
     Destructor.
   */
  ~PROTO_STEADY_DIFFUSION ();
  
  /**
     Assemble the conductivity matrix K(phi) (see the class comment).
  */
  bool assemble_conductivity_matrix ();
  /**
     Assemble the source (rhs) terms (see the class comment).
  */
  bool assemble_source_terms ();
  /**
     Solve for the unknown field phi.  (Both the conductivity matrix
     and the source terms must have been assembled prior to this call.)
  */
  bool solve ();
  /**
     Get the LES (Linear Equation Solver) that is being
     used by the protocol to solve the diffusion equation.
  */
  STEADY_DIFFUSION::LES *les () { return _les; }
  /**
     Get the primary field that this protocol is used
     to solve for.
  */
  FIELD_SCALAR *phi () { return _phi; }
  /**
     Get the geometry field that this protocol is supposed
     to use in the solution.
  */
  FIELD_VECTOR *geometry () { return _geometry; }
  /**
     Get the associated database.
  */
  DB *db () { return _algo->mgr()->db(); }
  /**
     Material.
  */
  MAT_DIFFUSION_ISO *mat (string mattype, string matname) ;
  
 private: // object data /////////////////////////////////////////////

  ALGO                                             *_algo;
  STEADY_DIFFUSION::LES                            *_les;
  FIELD_VECTOR                                     *_geometry;
  list <ECELL_STEADY_DIFFUSION *>                   _ecells;
  FIELD_SCALAR                                     *_phi;
  MAT_MAP <MAT_DIFFUSION_ISO>                      *_mat_map;

};

#endif
