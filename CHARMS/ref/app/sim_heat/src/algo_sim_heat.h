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
#ifndef ALGO_SIM_HEAT_H
# define ALGO_SIM_HEAT_H

#include <list>
#include <iterator>
#include "field.h"
#include "gcell.h"
#include "../../steady_diffusion/src/proto_steady_diffusion.h"
#include "les.h"
#include "algo.h"
#include "db.h"
#include "ebc.h"
#include "algo_refine.h"
#include "mgr.h"

class ALGO_SIM_HEAT : public ALGO{

 public: // object functions ////////////////////////////////////////

  ALGO_SIM_HEAT (string name, MGR *mgr, double tstart, double tend, double dt);

  ~ALGO_SIM_HEAT ();

  void setup ();
  
  void solve ();
  
  FIELD_VECTOR * geometry () {return _algo_steady_diffusion->geometry();} 
  
  FIELD_SCALAR *phi () { return _algo_steady_diffusion->phi(); }

  GMESH *gmesh () { return _algo_steady_diffusion->gmesh(); }

  void adapt ();

  private: // object data //////////////////////////////////////////
  
  MGR                    *_mgr; // reference only
  double                  _tstart, _tend, _dt, _time;
  ALGO_STEADY_DIFFUSION  *_algo_steady_diffusion ;
  ALGO_REFINE            *_algo_refine;
};

#endif
