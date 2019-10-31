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
#ifndef ALGO_ERR_REP_H
#define ALGO_ERR_REP_H

#include "fen.h"
#include "algo.h"
#include "gmesh.h"
#include "db.h"
#include "field.h"
#include "mgr.h"

/**
   This class defines a trivial interface for error reporting.
*/
class ALGO_ERR_REP : public ALGO {
  
 public: // object methods //@@//@@//@@//@@//@@//@@//@@//@@//@@//@@//
  
  ALGO_ERR_REP (string name, MGR *mgr) : ALGO (name, mgr) {}
  virtual ~ALGO_ERR_REP () {}
  /**
     Get `error' for given basis function.
   */
  virtual double errval (BFUN *bfun) { return 0; }
  /**
     This is for debugging.
  */
  virtual void write_mesh_sizeto_name (string filename) {}

};

#endif
