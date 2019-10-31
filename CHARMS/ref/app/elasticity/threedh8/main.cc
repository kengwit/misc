/***************************************************************************
                          main.cc  -  description
                             -------------------
    begin                : Sun Apr 1 2001
    copyright            : (C) 2001 by Petr Krysl
    email                : pkrysl@ucsd.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "stdlib.h"
#include "stdio.h"
#include "args.h"
#include "gmesh.h"
#include "gsubmesh.h"
#include "gcell_solid_h8.h"
#include "gcell_surf_q4.h"
#include "ecell_elasticity_h8.h"
#include "ecell_elasticity_h8_bbar.h"
#include "ecell_errest_elasticity_h8.h"
#include "ecell_elasticity_traction_q4.h"
#include "ecell_hexvue_h8.h"
#include "mat_elasticity_iso.h"
#include "load_on_gcell_body.h"
#include "mesh_mgr.h"
#include "algo_elasticity.h"
#include "algo_hexvue.h"
#include "algo_refine.h"
#include "algo_show_hier.h"

int
main (int argc, char **argv)
{
  ARGS::start_process_options (argc, argv, "algorithm name");
  int nadapt = ARGS::get_opt_int ("-nadapt", 3, "# number of adaptation cycles");
  int whier = ARGS::get_opt_int ("-whier", 0, "# write mesh hierarchy; true or false");
  int wdispmag = ARGS::get_opt_int ("-wdispmag", 0, "# write displacements + disp magnitude; true or false");
  int whvu = ARGS::get_opt_int ("-whvu", 1, "# write HEXVUE file; true or false");
  double node_size_mult = ARGS::get_opt_double ("-node_size_mult", 0.06, "node size multiplier");
  ARGS::finish_process_options ();
  
  if (argc < 2) {
    cerr << "Expected algorithm name as the first argument" << endl;
    exit (1);
  }
  string algorithm = argv[1];

  fles_initialize (&argc, &argv);

  GCELL_SOLID_H8::register_make_func ();
  GCELL_SURF_Q4::register_make_func ();
  ECELL_ELASTICITY_H8::register_make_func ();
  ECELL_ELASTICITY_H8_BBAR::register_make_func ();
  ECELL_ELASTICITY_TRACTION_Q4::register_make_func ();
  ECELL_HEXVUE_H8<3>::register_make_func ();
  ECELL_HEXVUE_H8<4>::register_make_func ();
  ECELL_ERREST_ELASTICITY_H8::register_make_func ();
  MAT_ELASTICITY_ISO::register_make_func ();
  LOAD_ON_GCELL_BODY::register_make_func();
  LOAD_ON_GCELL_TRACTION::register_make_func();

  DB db;
  db.load_file (algorithm + ".fpar");

  MGR mgr (&db);
  SMART_HANDLE<LOGGER_STREAM> toplsb = mgr.logger ("FAMULS (C) 2001-2002 P. Krysl", true);

  ALGO_ELASTICITY a (algorithm, &mgr);
  a.setup ();
  a.solve ();
  if (whvu) {
    if (wdispmag) {
      FIXED_VECTOR<4> zero(0);
      FIELD<4> *dispmag = new FIELD<4> ("dispmag", a.u()->bfun_set(), zero);
      for (int j = 0, n = dispmag->npairs(); j < n; j++) {
        FIELD_PAIR<4> *fp = dispmag->jth_field_pair (j);
        FIXED_VECTOR<4> v;
        FIXED_VECTOR<3> dp = a.u()->jth_field_pair(j)->dofparam();
        v(0) = dp(0);
        v(1) = dp(1);
        v(2) = dp(2);
        v(3) = dp.l2_norm();
        fp->set_dofparam (v);
      }
      ALGO_HEXVUE<4> h ("hexvue", &mgr);
      h.write_field (a.geometry(), a.geometry(), dispmag);
      delete dispmag;
    } else { // Only the displacements
      ALGO_HEXVUE<3> h ("hexvue", &mgr);
      h.write_field (a.geometry (), a.geometry (), a.u ());
    }
  }

  string progname = string (argv[0]);
  bool adaptive = false;
  string top = "algorithms/elasticity/" + algorithm;
  string path = top + "/adapt/adaptive";
  if (db.param_defined (path)) {
    adaptive = db.DB_GET_BOOL (path);
  }
  
  if ((progname == "adapt") || (adaptive)) {
    
  //  ALGO_REFINE r (a.name(), &mgr, a.gmesh());
    ALGO_SHOW_HIER sh (a.name(), &mgr);
    sh.set_node_size_mult (node_size_mult);

    for (int cycle = 0; cycle < nadapt; cycle++) {
      char buf[64];
      
      SMART_HANDLE<LOGGER_STREAM> lsb = mgr.logger("Adaptation cycle", true);
      *lsb << "cycle " << cycle << endl;
      
      a.adapt ();

      if (whier) {
        sprintf (buf, "-hier%d.egf", cycle);
        sh.save_as_elixir_file (algorithm + string (buf), a.u());
      }

      a.solve ();

      if (whvu) {
        if (wdispmag) {
          FIXED_VECTOR<4> zero(0);
          FIELD<4> *dispmag = new FIELD<4> ("dispmag", a.u()->bfun_set(), zero);
          for (int j = 0, n = dispmag->npairs(); j < n; j++) {
            FIELD_PAIR<4> *fp = dispmag->jth_field_pair (j);
            FIXED_VECTOR<4> v;
            FIXED_VECTOR<3> dp = a.u()->jth_field_pair(j)->dofparam();
            v(0) = dp(0);
            v(1) = dp(1);
            v(2) = dp(2);
            v(3) = dp.l2_norm();
            fp->set_dofparam (v);
          }
          sprintf (buf, "-adapt%d.hvu", cycle);
          ALGO_HEXVUE<4> h ("hexvue", &mgr);
          h.write_field_to_name (algorithm + string (buf),
                                 a.geometry(), a.geometry(), dispmag);
          delete dispmag;
        } else { // Only the displacements
          ALGO_HEXVUE<3> h ("hexvue", &mgr);
          sprintf (buf, "-adapt%d.hvu", cycle);
          h.write_field_to_name (algorithm + string (buf),
                                 a.geometry(), a.geometry(), a.u());
        }
      }
      
      *lsb << "done in adaptation cycle " << cycle << endl;
      
    }
    
  }

  fles_finalize ();

  *toplsb << "done" << endl;

  return 0;
}
