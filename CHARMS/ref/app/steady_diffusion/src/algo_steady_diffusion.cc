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
#include "algo_steady_diffusion.h"
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

ALGO_STEADY_DIFFUSION::ALGO_STEADY_DIFFUSION (string name, MGR *mgr)
  : ALGO::ALGO (name, mgr)
{
  _gmesh = 0;
  CHECK (mgr != 0, EXCEPTION_NULL_PTR,;);
  string top = "algorithms/steady_diffusion/" + this->name ();
  string gmesh_name = this->mgr()->db()->DB_GET_STRING (top + "/gmesh");
  _gmesh = this->mgr()->mesh_mgr()->gmesh (gmesh_name);
  _gsubmeshes.clear();
  _algo_refine = new ALGO_REFINE(name, mgr, _gmesh);
  _geometry = 0;
  _phi = 0;
  _ebc = 0;
  _shape_approx = 0;
  _algo_errest = 0;
  _time = 0;
  
}

void
ALGO_STEADY_DIFFUSION::setup ()
{
  CHECK (this->mgr()->db() != 0, EXCEPTION_NULL_PTR,;);
  string top = "algorithms/steady_diffusion/" + name ();

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
  // SMART_HANDLE<BFUN_SET> geom_bfun_set (BFUN_SET::make (_gmesh, 0)); //RAW
  SMART_HANDLE<BFUN_SET> geom_bfun_set (BFUN_SET::make (_algo_refine->refinement_strategy(),_gmesh)); //RAW
  _geometry = new FIELD_VECTOR ("geometry", geom_bfun_set, zero);
  ALGO_GEOM init_geom ("", this->mgr());
  init_geom.init (geometry(), _gmesh);

  // create the phi field
  // SMART_HANDLE<BFUN_SET> phi_bfun_set (BFUN_SET::make (_gsubmeshes, 0)); //RAW
  SMART_HANDLE<BFUN_SET> phi_bfun_set (BFUN_SET::make (_algo_refine->refinement_strategy(),_gsubmeshes)); //RAW
  _phi      = new FIELD_SCALAR ("phi", phi_bfun_set, 0);

  // read the essential boundary conditions
  string path = top + "/essential_bcs/ebc_file";
  string file = this->mgr()->db()->DB_GET_STRING (path);
  _ebc = new EBC<1> (_gmesh);
  CHECK (_ebc->read (file), EXCEPTION_BAD_VALUE,;);
  if (this->mgr()->db()->param_defined (top + "/shape_approx/shp_file")) {
    string path = top + "/shape_approx/shp_file";
    string file = this->mgr()->db()->DB_GET_STRING (path);
    _shape_approx = new SHAPE_APPROX<3> (_gmesh);
    CHECK (_shape_approx->read (file), EXCEPTION_BAD_VALUE,;);
  }

}

void
ALGO_STEADY_DIFFUSION::solve ()
{
  CHECK (this->mgr()->db() != 0, EXCEPTION_NULL_PTR,;);
  string top = "algorithms/steady_diffusion/" + name ();

  SMART_HANDLE<LOGGER_STREAM> lsb = this->mgr()->logger("solve", true);
  
  // Apply EBC to the unknown field
  *lsb << "applying ebc ... ";
  _ebc->apply (this->phi());
   *lsb << "done" << endl;

  if (_shape_approx ) {
    *lsb << "applying shape approx ... ";
    _shape_approx->apply (this->geometry());
    *lsb << "done" << endl;
   }

  // Set up solver and protocol
  string les_name = this->mgr()->db()->DB_GET_STRING (top + "/les");
  bool use_nnzs = true;
  string use_nnzs_param = "les/" + les_name + "/use_nnzs";
  if (this->mgr()->db()->param_defined (use_nnzs_param)) {
    use_nnzs = this->mgr()->db()->DB_GET_BOOL (use_nnzs_param);
  }
  size_t max_per_row = 500;
  string max_per_row_param = "les/" + les_name + "/max_per_row";
  if (this->mgr()->db()->param_defined (max_per_row_param)) {
    max_per_row = this->mgr()->db()->DB_GET_INTEGER (max_per_row_param);
  }
  double nz_per_row = (double) max_per_row / (double) this->phi()->ndofparams();
  *lsb << "computing " << this->phi()->ndofparams() << " degrees of freedom" << endl;
  if (nz_per_row > 1.0) nz_per_row = 1.0;
  char options[512];
  if (use_nnzs) {
    *lsb << "using per-row non-zeros" << endl;
    sprintf (options, "use_nnzs");
  } else {
    *lsb << "estimated nz_per_row as " << nz_per_row <<
      " (" << max_per_row << " / " << this->phi()->ndofparams() << ")" << endl;
    sprintf (options, "nz_per_row %g", nz_per_row);
  }
  STEADY_DIFFUSION::LES les (les_name, this->mgr(), string (options));
  PROTO_STEADY_DIFFUSION proto_steady_diffusion (this, &les, geometry(), this->phi());
  // Construct conductivity matrix and the source term (rhs)
  *lsb << "assembling conductivity matrix ... ";
  proto_steady_diffusion.assemble_conductivity_matrix ();
  *lsb << "done" << endl;
  *lsb << "assembling source term ... ";
  proto_steady_diffusion.assemble_source_terms ();
  *lsb << "done" << endl;
  *lsb << "solving the linear system ... ";
  proto_steady_diffusion.solve ();
  *lsb << "done" << endl;
}

 
void
ALGO_STEADY_DIFFUSION::adapt ()
{
  const string phi_name = "phi";
  const string geometry_name = "geom";
  SMART_HANDLE<LOGGER_STREAM> lsb = this->mgr()->logger("adapt", true);

  if (_algo_errest == 0) {
    string path = "algorithms/refine/" + this->name();
    string errestim = "gradation";
    if (this->mgr()->db()->param_defined (path + "/errestim")) {
      errestim = this->mgr()->db()->DB_GET_STRING (path + "/errestim");
    }
    if (errestim == "gradation") {
      _algo_errest = new ALGO_ERREST_GRADATION<1> (this->name(), this->mgr(), _gmesh);
    } else {
      _algo_errest = new ALGO_ERREST_ENERGY<1> (this->name(), this->mgr());
    }
  }
  
  *lsb << "assembling error ... ";
  // Assemble the error, and then use it to update the refinement context
  _algo_errest->assemble_error (geometry(), this->phi(), _time) ;
 
  *lsb << "done" << endl;
  *lsb << "updating refinement context ... ";
  SMART_HANDLE <BFUN_SET> new_bfun_set = _algo_refine->adapt (_algo_errest, this->phi());
  //new_bfun_set->debug_display("junk.jpp");//RAW
  *lsb << "done" << endl;

  *lsb << "making new fields ... ";
  FIELD_SCALAR *new_phi = new FIELD_SCALAR (phi_name, new_bfun_set, 0);
  

  FIXED_VECTOR<3>  zero(0);
  FIELD_VECTOR *new_geometry = new FIELD_VECTOR (geometry_name, new_bfun_set, zero);
  
   
  *lsb << "done" << endl;
    
  *lsb << "transfering phi ... ";
  PROTO_FIELD_TRANSFER<1> transf1;
  transf1.transfer (this->phi(), new_phi);
   *lsb << "done" << endl;
  *lsb << "transfering geometry ... ";
  PROTO_FIELD_TRANSFER<3> transf3;
  transf3.transfer (geometry(), new_geometry);
  *lsb << "done" << endl;
  
  {
    //BFUN_SET *bfun_set = _phi->bfun_set();
    delete _phi;
    // delete bfun_set; RAW the bfun_set should not be deleted by calling its destructor, but rather released
    // Maybe a handle would be a good solution
  }
  _phi = new_phi;
  //_phi->bfun_set()->debug_display("junk_phi.jpp");//RAW
   {
    //BFUN_SET *bfun_set = _geometry->bfun_set();
    delete _geometry;
    // delete bfun_set; RAW see discussion above
  }
  _geometry = new_geometry;
  //_geometry->bfun_set()->debug_display("junk_geometry.jpp");//RAW
   *lsb << "done" << endl;
}

ALGO_STEADY_DIFFUSION::~ALGO_STEADY_DIFFUSION () {
  if (_geometry) delete _geometry;
  if (_phi) delete _phi;
  if (_algo_errest) delete _algo_errest;
  if (_ebc) delete _ebc;
}

void
ALGO_STEADY_DIFFUSION::adapt (double time)
{
   _time = time;
  adapt();
  _time = 0;
}
