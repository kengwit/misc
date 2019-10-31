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
#include "ecell_wave_h8_opt.h"
// #include "ecell_wave_h8_damp.h"
#include "ecell_wave_traction_q4.h"
#include "ecell_errest_wave_h8.h"
#include "ecell_hexvue_h8.h"
#include "ecell_xhexvue_h8.h"
#include "mat_viscel_iso.h"
#include "load_on_gcell_body.h"
#include "mesh_mgr.h"
#include "algo_wave.h"
#include "algo_hexvue.h"
#include "algo_refine.h"
#include "algo_show_hier.h"
extern "C" {
#include "fles.h"
}


int
main (int argc, char **argv)
{
  ARGS::start_process_options (argc, argv, "algorithm name");
  int nadapt = ARGS::get_opt_int ("-nadapt", 3, "# number of adaptation cycles");
  int whier = ARGS::get_opt_int ("-whier", 0, "# write mesh hierarchy; true or false");
  int wdispmag = ARGS::get_opt_int ("-wdispmag", 0, "# write displacements + disp magnitude; true or false");
  int whvu = ARGS::get_opt_int ("-whvu", 1, "# write HEXVUE file; true or false");
  double node_size_mult = ARGS::get_opt_double ("-node_size_mult", 0.06, "node size multiplier");
  string hist_file_prefix = ARGS::get_opt_string ("-hist_file_prefix", "",
                                                 "# history file prefix");
  ARGS::finish_process_options ();
  
  if (argc < 2) {
    cerr << "Expected algorithm name as the first argument" << endl;
    exit (1);
  }
  string algorithm = argv[1];

  fles_initialize (&argc, &argv); //RAW

  GCELL_SOLID_H8::register_make_func ();
  ECELL_WAVE_H8_OPT::register_make_func ();
//  ECELL_WAVE_H8_DAMP::register_make_func ();
  GCELL_SURF_Q4::register_make_func ();
  ECELL_WAVE_TRACTION_Q4::register_make_func ();
  ECELL_HEXVUE_H8<1>::register_make_func ();
  ECELL_HEXVUE_H8<3>::register_make_func ();
  ECELL_HEXVUE_H8<4>::register_make_func ();
  ECELL_XHEXVUE_H8::register_make_func ();
  ECELL_ERREST_WAVE_H8::register_make_func ();
  MAT_VISCEL_ISO::register_make_func ();
  LOAD_ON_GCELL_BODY::register_make_func();
  LOAD_ON_GCELL_TRACTION::register_make_func();

  DB db;
  db.load_file (algorithm + ".fpar");

  MGR mgr (&db);
  SMART_HANDLE<LOGGER_STREAM> toplsb = mgr.logger ("FAMULS (C) 2001-2005 P. Krysl", true);

  ALGO_WAVE a (algorithm, &mgr);
  a.set_hist_file_prefix (hist_file_prefix);
  a.setup ();

  ALGO_SHOW_HIER sh (a.name(), &mgr);
  sh.set_node_size_mult (node_size_mult);
  
  a.set_nadapt (nadapt);
  a.run ();
  
  *toplsb << "done" << endl;

  fles_finalize ();

  return 0;
}
