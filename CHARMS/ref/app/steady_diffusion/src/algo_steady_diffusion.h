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
#ifndef ALGO_STEADY_DIFFUSION_H
# define ALGO_STEADY_DIFFUSION_H

#include <list>
#include <iterator>
#include "field.h"
#include "gcell.h"
#include "proto_steady_diffusion.h"
#include "select_les.h"
#include "algo.h"
#include "db.h"
#include "ebc.h"
#include "shape_approx.h"
#include "algo_refine.h"
#include "mgr.h"

/**
   This class implements the steady diffusion algorithm.
   It is assumed that the embedding space dimension is equal to the manifold
   dimension of the cells.  In other words, the steady diffusion equation is
   written in the Euclidean space R^n.
*/
class ALGO_STEADY_DIFFUSION : public ALGO {

 public: // object functions ////////////////////////////////////////

  /**
     Construct an algorithm.  Give it a name, and let it load
     its parameters from the database.  The mesh manager is used
     to create the gmesh specified for the algorithm in the database.
  */
  ALGO_STEADY_DIFFUSION (string name, MGR *mgr);
  
  /**
   */
  ~ALGO_STEADY_DIFFUSION ();
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
  FIELD_VECTOR * geometry () {return _geometry;}
  /**
     Return the primary field.
  */
  FIELD_SCALAR *phi () { return _phi; }
  /**
     Return the mesh.
  */
  GMESH *gmesh () { return _gmesh; }
  /**
     Adapt the basis function sets for the primary field and the geometry.
  */
  void adapt ();

/**
     Same as adapt but increase in internal time takes place each time when function is being called.  
*/
  void adapt (double time);
  

 private: // object functions /////////////////////////////////////

 private: // object data //////////////////////////////////////////

  MGR                    *_mgr; // reference only
  GMESH                  *_gmesh; // reference only
  list <GSUBMESH *>       _gsubmeshes; // references only
  FIELD_VECTOR           *_geometry; // responsible for 
  FIELD_SCALAR           *_phi; // responsible for
  EBC<1>                 *_ebc; // responsible for
  SHAPE_APPROX<3>        *_shape_approx; // responsible for
  ALGO_ERREST<1>         *_algo_errest; // responsible for
  ALGO_REFINE            *_algo_refine; 
  double                 _time;    
};

#endif
