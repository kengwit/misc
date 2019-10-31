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
#include "fen_maker.h"
#include "gmesh.h"
#include "mesh_mgr.h"
#include "gsubmesh.h"
#include "ref_ctx.h"
#include "bfun_set.h"
#include "../../steady_diffusion/src/algo_steady_diffusion.h"
#include "fixed_vector.h"
#include "algo_geom.h"
#include "ebc.h"
#include "algo_refine.h"
#include "proto_field_transfer.h"
#include "algo_errest_gradation.h"
#include "algo_errest_energy.h"
#include "algo_xmgr.h"
#include "logger_stream.h"
#include "algo_hexvue.h"
#include "algo_sim_heat.h"

ALGO_SIM_HEAT::ALGO_SIM_HEAT (string name, MGR *mgr, double tstart, double tend, double dt): ALGO::ALGO (name, mgr)
{
  _tstart = tstart;
  _tend   = tend;
  _dt     = dt;
  _time   = _tstart; 
  _algo_steady_diffusion = new ALGO_STEADY_DIFFUSION (name, mgr);
  
}

void ALGO_SIM_HEAT::setup () {_algo_steady_diffusion->setup();}

void ALGO_SIM_HEAT::solve () { _algo_steady_diffusion->solve();}

 
void
ALGO_SIM_HEAT::adapt ()
{
  _time = _time + _dt;
  _algo_steady_diffusion->adapt(_time);   
}

ALGO_SIM_HEAT::~ALGO_SIM_HEAT () {}
  
