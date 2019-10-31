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
#ifndef ALGO_MODL_H
# define ALGO_MODL_H

#include <list>
#include <iterator>
#include "field.h"
#include "gcell.h"
#include "proto_modl.h"
#include "evs_ooofs.h"
#include "algo.h"
#include "db.h"
#include "ebc.h"
#include "algo_refine.h"
#include "mesh_mgr.h"

class ALGO_MODL : public ALGO {

 public: // object functions ////////////////////////////////////////

  /**
     Construct an algorithm.  Give it a name, and let it load
     its parameters from the database.  The mesh given on input
     must be the one specified for the algorithm in the database: an
     exception is raised if that is not the case.
  */
  ALGO_MODL (string name, MGR *mgr);
  /**
   */
  ~ALGO_MODL ();
  /**
     Prepare the algorithm for solution.  
  */
  void setup ();
  /**
     Solve for the primary field.
  */
  void solve ();
  /**
     Return the geometry.
  */
  FIELD_VECTOR * geometry () { return _geometry; }
  /**
     Return the primary field (displacements).
  */
  FIELD_VECTOR *u () { return _u; }
  /**
     Return the mesh.
  */
  GMESH *gmesh () { return _gmesh; }
  /**
     Adapt the basis function sets for the primary field and the geometry.
  */
  void adapt ();



 private: // object data //////////////////////////////////////////

  GMESH                  *_gmesh;
  list <GSUBMESH *>       _gsubmeshes;
  FIELD_VECTOR           *_geometry;
  FIELD_VECTOR           *_u;
  EBC<3>                 *_ebc;
  ALGO_ERREST<3>         *_algo_errest;
  ALGO_REFINE            *_algo_refine;
  double                  _shift;
  bool                    _is_full_stress_tensor;
  void write_modes (EVS_OOOFS &evs, PROTO_MODL &proto_modl);
  void report_strain_energy_per_gcell_group (PROTO_MODL &proto_modl);

};

#endif




