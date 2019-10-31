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
#include "logger_stream_mgr.h"
#include "logger_stream_op.h"
#include "bfun_set.h"
#include "algo_gradation.h"
#include "algo_geom.h"
#include "algo_hexvue.h"
#include "field.h"
#include "gcell.h"

ALGO_GRADATION::ALGO_GRADATION (string name, MGR *mgr, GMESH *gmesh,double t) : ALGO (name, mgr)
{
  string path = "algorithms/gradation/" + this->name();
  _gmesh = gmesh;
  _mesh_size = 0;
  _geometry = 0;
  
  string expr[1];
  expr[0] = this->mgr()->db()->DB_GET_STRING (path + "/mesh_size_func");
  _h = new FUNC<1> (expr);
  // cerr << "Mesh size func " << expr[0] << endl;

  _bfun_set = BFUN_SET::make (BFUN_SET::QUASI_HIERARCHICAL,gmesh); 
  _mesh_size = new FIELD_SCALAR (this->name(), _bfun_set, 0.0); 
  init_mesh_size(t);
  POINT zero (0.0);
  _geometry = new FIELD_VECTOR (this->name(), _bfun_set, zero); 
  ALGO_GEOM init_geom ("", this->mgr());
  init_geom.init (_geometry, gmesh);
}

void ALGO_GRADATION::init_mesh_size(double t)
{
  POINT at1; POINT4 at4;
  FIXED_VECTOR<1> hval;
  for (sizet j = 0; j < _mesh_size->npairs (); j++) {
    FIELD_PAIR<1> *fp = _mesh_size->field_pair (j);
    FEN *fen = fp->bfun()->fen();
    at1 = fen->ref_loc();
    at4(0) = at1(0); at4(1) = at1(1);at4(2) = at1(2); at4(3) = t;
    hval = (*_h)(at4);
    fp->set_dofparam (hval);
  }
}


ALGO_GRADATION::~ALGO_GRADATION ()
{
  if (_mesh_size) delete _mesh_size;
  if (_geometry) delete _geometry;
  if (_h) delete _h;
}


double 
ALGO_GRADATION::mesh_size (EVALPT &evalpt)
{
  GCELL *gcell = evalpt.gcell ();
  CHECK(gcell != 0,EXCEPTION_NULL_PTR,;);
  POINT mapped_to_param_loc;
  POINT param_loc = evalpt.param_loc();
  GCELL *evgcell
    = gcell->map_to_topmost_gcell_in_field (_mesh_size, param_loc,  &mapped_to_param_loc);
  CHECK(evgcell != 0,EXCEPTION_NULL_PTR,;);
  FIXED_VECTOR<1> result = _mesh_size->evaluate (_geometry, evgcell, mapped_to_param_loc);
  return result(0);
}


void
ALGO_GRADATION::write_mesh_size_to_name (string filename)
{
  ALGO_HEXVUE<1> h("hexvue", this->mgr());
  h.write_field_to_name (filename, _geometry, _geometry, _mesh_size);
}
  
