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
#ifndef ALGO_GEOM_H
#define ALGO_GEOM_H

#include "algo.h"
#include "fen.h"
#include "gmesh.h"
#include "mgr.h"

/**
  *@author Petr Krysl
  */

class ALGO_GEOM : public ALGO {
  
 public: 
  
  /**
     This algorithm works on geometry fields.
   */
  ALGO_GEOM (string name, MGR *mgr);
  ~ALGO_GEOM();
  /**
     Initialize the geometry field from the reference geometry of the mesh.
     This operation may be performed *only* on the initial mesh, not when
     the mesh had been refined.  Use the interpolation algorithm for fields
     on refined meshes.  
  */
  void init (FIELD_VECTOR *geometry, GMESH *gmesh);
  
};

#endif
