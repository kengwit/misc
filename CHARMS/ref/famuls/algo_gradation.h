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
#ifndef ALGO_GRADATION_H
#define ALGO_GRADATION_H

#include "fen.h"
#include "algo.h"
#include "gmesh.h"
#include "db.h"
#include "field.h"
#include "func.h"
#include "evalpt.h"
#include "smart_handle.h"

class ALGO_GRADATION : public ALGO {
  
 public: // object methods //@@//@@//@@//@@//@@//@@//@@//@@//@@//@@//
  
  ALGO_GRADATION (string name, MGR *mgr, GMESH *gmesh, double t);
  ~ALGO_GRADATION ();

  /**
    This has to be called to initilize mesh size field for a time dependent mesh 
    size function.
  */ 
  void init_mesh_size(double t);
  double mesh_size (EVALPT &evalpt);
  void write_mesh_size_to_name (string filename);

 private: // object data //@@//@@//@@//@@//@@//@@//@@//@@//@@//@@//@@

  GMESH                  *_gmesh;
  SMART_HANDLE<BFUN_SET>  _bfun_set;
  FIELD<1>               *_mesh_size;
  FIELD_VECTOR           *_geometry;
  FUNC<1>                *_h;

};

#endif
