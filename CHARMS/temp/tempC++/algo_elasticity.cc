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
#include "gsubmesh.h"
#include "ref_ctx.h"
#include "bfun_set.h"
#include "algo_elasticity.h"
#include "fixed_vector.h"
#include "algo_geom.h"
#include "ebc.h"
#include "algo_errest_gradation.h"
#include "algo_errest_energy.h"
#include "algo_refine.h"
#include "proto_field_transfer.h"
#include "algo_xmgr.h"
#include "logger_stream.h"
#include "smart_handle.h"

ALGO_ELASTICITY::ALGO_ELASTICITY (string name, MGR *mgr) : ALGO (name, mgr)
{
  CHECK (this->mgr()->db() != 0, EXCEPTION_NULL_PTR,;);
  _gmesh = 0;
  string top = "algorithms/elasticity/" + this->name ();
  string gmesh_name = this->mgr()->db()->DB_GET_STRING (top + "/gmesh");
  _gmesh = this->mgr()->mesh_mgr()->gmesh (gmesh_name);
  _gsubmeshes.clear();
  _algo_refine = new ALGO_REFINE(name, mgr, _gmesh);
  _geometry = 0;
  _u = 0;
  _ebc = 0;
  _algo_errest = 0;
}

void
ALGO_ELASTICITY::setup ()
{
  CHECK (this->mgr()->db() != 0, EXCEPTION_NULL_PTR,;);
  string top = "algorithms/elasticity/" + name ();

  // construct the gsubmesh list
  list <string> gsubmesh_names = this->mgr()->db()->DB_GET_STRING_LIST (top + "/gsubmeshes");
  {
    list <string>::const_iterator i = gsubmesh_names.begin ();
    while (i != gsubmesh_names.end ()) {
      string n = *i;
      GSUBMESH *gsubmesh = _gmesh->find_gsubmesh (n);
      CHECK (gsubmesh, EXCEPTION_NULL_PTR,;);
      _gsubmeshes.push_back (gsubmesh);
      i++;
    }
  }

  // create the geometry field
  FIXED_VECTOR<3>  zero(0);
  SMART_HANDLE<BFUN_SET> geom_bfun_set = BFUN_SET::make (_algo_refine->refinement_strategy(),_gmesh);
  _geometry = new FIELD_VECTOR ("geometry", geom_bfun_set, zero);
  ALGO_GEOM init_geom ("", this->mgr());
  init_geom.init (_geometry, _gmesh);

  // create the field
  SMART_HANDLE<BFUN_SET> u_bfun_set = BFUN_SET::make (_algo_refine->refinement_strategy(),_gsubmeshes);
  _u = new FIELD_VECTOR ("u", u_bfun_set, zero);

  // read the essential boundary conditions
  string path = top + "/essential_bcs/ebc_file";
  string file = this->mgr()->db()->DB_GET_STRING (path);
  _ebc = new EBC<3> (_gmesh);
  CHECK (_ebc->read (file), EXCEPTION_BAD_VALUE,;);
}

void
ALGO_ELASTICITY::solve ()
{
  CHECK (this->mgr()->db() != 0, EXCEPTION_NULL_PTR,;);
  string top = "algorithms/elasticity/" + name ();

  SMART_HANDLE<LOGGER_STREAM> lsb = this->mgr()->logger("solve", true);

  // Apply EBC to the unknown field
  *lsb << "applying EBC's ... ";
  _ebc->apply(_u);
  //_u->attach_ebc(_ebc);
  //_u->enforce_ebc();
  
  *lsb << "done" << endl;
  
  // Set up solver and protocol
  string les_name = this->mgr()->db()->DB_GET_STRING (top + "/les");
  bool use_nnzs = true;
  string use_nnzs_param = "les/" + les_name + "/use_nnzs";
  if (this->mgr()->db()->param_defined (use_nnzs_param)) {
    use_nnzs = this->mgr()->db()->DB_GET_BOOL (use_nnzs_param);
  }
  size_t max_per_row = 1000;
  string max_per_row_param = "les/" + les_name + "/max_per_row";
  if (this->mgr()->db()->param_defined (max_per_row_param)) {
    max_per_row = this->mgr()->db()->DB_GET_INTEGER (max_per_row_param);
  }
  double nz_per_row = (double) max_per_row / (double) _u->ndofparams();
  *lsb << "computing " << _u->ndofparams() << " degrees of freedom" << endl;
  if (nz_per_row > 1.0) nz_per_row = 1.0;
  char options[512];
  if (use_nnzs) {
    *lsb << "using per-row non-zeros" << endl;
    sprintf (options, "use_nnzs");
  } else {
    *lsb << "estimated nz_per_row as " << nz_per_row <<
      " (" << max_per_row << " / " << _u->ndofparams() << ")" << endl;
    sprintf (options, "nz_per_row %g", nz_per_row);
  }
  LES les (les_name, this->mgr(), string (options));
  PROTO_ELASTICITY proto_elasticity (this, this->mgr()->db(), &les, _geometry, _u);

  // Now assemble the lhs, rhs, and solve
  *lsb << "assembling stiffness matrix ... ";
  proto_elasticity.assemble_stiffness_matrix ();
  *lsb << "done"  << endl;
  *lsb << "assembling load terms ... ";
  proto_elasticity.assemble_load_terms ();
  *lsb << "done"  << endl;
  *lsb << "solving the linear system ... ";
  proto_elasticity.solve ();
  *lsb << "done"  << endl;
}

 
void
ALGO_ELASTICITY::adapt ()
{
  string u_name = "u";
  string geometry_name = "geom";

  SMART_HANDLE<LOGGER_STREAM> lsb = this->mgr()->logger("adapt", true);

  if (_algo_errest == 0) {
    string path = "algorithms/refine/" + this->name();
    string errestim = "gradation";
    if (this->mgr()->db()->param_defined (path + "/errestim")) {
      errestim = this->mgr()->db()->DB_GET_STRING (path + "/errestim");
    }
    if (errestim == "gradation") {
      _algo_errest = new ALGO_ERREST_GRADATION<3> (this->name(), this->mgr(), _gmesh);
    } else {
      _algo_errest = new ALGO_ERREST_ENERGY<3> (this->name(), this->mgr());
    }
  }
  
  *lsb << "assembling error ... ";
  // Assemble the error, and then use it to update the refinement context
  _algo_errest->assemble_error (_geometry, _u);
  *lsb << "done" << endl;
  *lsb << "updating refinement context ... ";
  SMART_HANDLE <BFUN_SET> new_bfun_set = _algo_refine->adapt(_algo_errest, _u); 
  *lsb << "done" << endl;

  *lsb << "making new fields ... ";
  FIXED_VECTOR<3>  zero(0);
  FIELD_VECTOR *new_u = new FIELD_VECTOR (u_name, new_bfun_set, zero);
  FIELD_VECTOR *new_geometry = new FIELD_VECTOR (geometry_name,new_bfun_set , zero);
  *lsb << "done" << endl;
    
  *lsb << "transfering u ... ";
  PROTO_FIELD_TRANSFER<3> transf3;
  transf3.transfer (_u, new_u);
  *lsb << "done" << endl;
  *lsb << "transfering geometry ... ";
  transf3.transfer (_geometry, new_geometry);
  *lsb << "done" << endl;
  
  delete _u;
  _u = new_u;
  delete _geometry;
  _geometry = new_geometry;
  
  //  delete new_u; //RAW
  //delete new_geometry; //RAW 
  *lsb << "done" << endl;
}

ALGO_ELASTICITY::~ALGO_ELASTICITY () {
  if (_geometry) delete _geometry;
  if (_u) delete _u;
  if (_algo_errest) delete _algo_errest;
  if (_ebc) delete _ebc;
}
