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
#include "algo_modl.h"
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
#include "ecell_hexvue_t10.h"
#include "algo_hexvue.h"


ALGO_MODL::ALGO_MODL (string name, MGR *mgr) : ALGO (name, mgr)
{
  CHECK (this->mgr()->db() != 0, EXCEPTION_NULL_PTR,;);
  _gmesh = 0;
  string top = "algorithms/modl/" + this->name ();
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
ALGO_MODL::setup ()
{
  CHECK (this->mgr()->db() != 0, EXCEPTION_NULL_PTR,;);
  string top = "algorithms/modl/" + name ();

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
ALGO_MODL::solve ()
{
  CHECK (this->mgr()->db() != 0, EXCEPTION_NULL_PTR,;);
  string top = "algorithms/modl/" + name ();

  SMART_HANDLE<LOGGER_STREAM> lsb = this->mgr()->logger("solve", true);

  // Apply EBC to the unknown field
  *lsb << "applying EBC's ... ";
  _ebc->apply(_u);
  //_u->attach_ebc(_ebc);
  //_u->enforce_ebc();
  
  *lsb << "done" << endl;

  // Set up solver and protocol
  string evs_name = this->mgr()->db()->DB_GET_STRING (top + "/evs");
  _shift = 0;
  string shift_s = "evs/" + evs_name + "/shift";
  if (this->mgr()->db()->param_defined (shift_s)) {
    _shift = this->mgr()->db()->DB_GET_DOUBLE (shift_s);
  }
  string fsn = "evs/" + evs_name + "/freq_shift";
  if (this->mgr()->db()->param_defined (fsn)) {
    double freq_shift = this->mgr()->db()->DB_GET_DOUBLE (fsn);
    _shift = pow(2*M_PI*freq_shift, 2);
  }
  int nmodes = 5;
  string nmodes_param = "evs/" + evs_name + "/nmodes";
  if (this->mgr()->db()->param_defined (nmodes_param)) {
    nmodes = this->mgr()->db()->DB_GET_INTEGER (nmodes_param);
  }
  double rtolv = 0.001;
  string rtolv_param = "evs/" + evs_name + "/rtolv";
  if (this->mgr()->db()->param_defined (rtolv_param)) {
    rtolv = this->mgr()->db()->DB_GET_DOUBLE (rtolv_param);
  }
  sizet max_iter = 40;
  string max_iter_param = "evs/" + evs_name + "/max_iter";
  if (this->mgr()->db()->param_defined (max_iter_param)) {
    max_iter = this->mgr()->db()->DB_GET_INTEGER (max_iter_param);
  }

  double init_p = 0;
  string pressure_param = "initial_cond/pressure";
  if (this->mgr()->db()->param_defined (pressure_param)) {
    init_p = this->mgr()->db()->DB_GET_DOUBLE (pressure_param);
  } 
  _is_full_stress_tensor = false;
  if (this->mgr()->db()->param_defined ("initial_cond/is_full_stress_tensor")) {
     _is_full_stress_tensor  = this->mgr()->db()->DB_GET_BOOL ("initial_cond/is_full_stress_tensor"); 
  }
	

  *lsb << "computing " << _u->ndofparams() << " degrees of freedom" << endl;
  char options[512];
  EVS_OOOFS evs(evs_name, this->mgr(), string (options));
  evs.set_nmodes(nmodes);
  evs.set_rtolv(rtolv);
  evs.set_max_iter(max_iter);
  PROTO_MODL proto_modl (this, this->mgr()->db(), &evs, _geometry, _u);

  // Now assemble the lhs, rhs, and solve
  *lsb << "applying shift " << _shift << endl;
  *lsb << "assembling stiffness matrix ... ";
  if ( !_is_full_stress_tensor) {
  proto_modl.assemble_stiffness_matrix (init_p, _shift);
  } else {
  proto_modl.assemble_stiffness_matrix ( _shift);
  } 
  *lsb << "done"  << endl;
  *lsb << "assembling mass matrix ...";
  proto_modl.assemble_mass_matrix ();
  *lsb << "total mass = " << proto_modl.total_mass() << endl;
  *lsb << "done"  << endl;
  *lsb << "solving the eigen value problem... ";
  proto_modl.solve (_shift);
  *lsb << "done"  << endl;
  *lsb << "writing mode shapes... ";
  write_modes (evs, proto_modl);
  *lsb << "done"  << endl;
}

void
ALGO_MODL::write_modes (EVS_OOOFS &evs, PROTO_MODL &proto_modl)
{
  FIXED_VECTOR<4> zero(0);
  SMART_HANDLE<LOGGER_STREAM> lsb = this->mgr()->logger("write modes", true);
  for (sizet mode = 0; mode < evs.nmodes(); mode++) {
    double eigval = evs.solution ("u", _u,  mode);
    double freq = sqrt(eigval+_shift)/2/M_PI;
    if (_shift == 0) { 
      *lsb << "Mode: " << mode+1 << "; frequency: " << freq << " [Hz] " << "; angular frequency: " << freq*2*M_PI << " [rad/s]" << endl;
    } else {
      *lsb << "Mode: " << mode+1 << "; frequency: " << freq << " (shifted: " << eigval << ")" << " [Hz] " << freq*2*M_PI << " [rad/s]" << endl;
    }
    report_strain_energy_per_gcell_group (proto_modl);
    char buf[1024];
    sprintf(buf, "mode:%d:freq:%g", mode+1, freq);
    FIELD<4> *dispmag = new FIELD<4> (buf, _u->bfun_set(), zero);
    for (int j = 0, n = dispmag->npairs(); j < n; j++) {
      FIELD_PAIR<4> *fp = dispmag->jth_field_pair (j);
      FIXED_VECTOR<4> v;
      FIXED_VECTOR<3> dp = _u->jth_field_pair(j)->dofparam();
      v(0) = dp(0);
      v(1) = dp(1);
      v(2) = dp(2);
      v(3) = dp.l2_norm();
      fp->set_dofparam (v);
    }
    ALGO_HEXVUE<4> h ("hexvue", this->mgr());
    sprintf (buf, "-mode%d.hvu", mode+1);
    h.write_field_to_name (this->name() + string (buf),_geometry, _geometry, dispmag);
    delete dispmag;
  }
  *lsb << "done"  << endl;
}


void
ALGO_MODL::report_strain_energy_per_gcell_group (PROTO_MODL &proto_modl)
{
  SMART_HANDLE<LOGGER_STREAM> lsb = this->mgr()->logger("write strain energies", true);
  GMESH::gsubmesh_enumerator_t gsme = _gmesh->gsubmesh_enumerator ();
  gsme.reset();
  GSUBMESH *gsm = 0;
  while (gsm = gsme.next()) {
    GSUBMESH::gcell_group_enumerator_t gge = gsm->gcell_group_enumerator ();
    gge.reset();
    GCELL_GROUP *gg = 0;
    while (gg = gge.next()) {
      *lsb << "Energy of gcell group " << gg->name() << ":\t"<<
        proto_modl.total_strain_energy_in_gcell_group(_u, gg) << endl;
    }
  }
  *lsb << "done"  << endl;
}



ALGO_MODL::~ALGO_MODL () {
  if (_geometry) delete _geometry;
  if (_u) delete _u;
  if (_algo_errest) delete _algo_errest;
  if (_ebc) delete _ebc;
}
