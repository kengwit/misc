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
#include "ecell_modl_h8.h"
//#include "ecell_modl_h8_bbar.h"
#include "ecell_hexvue_h8.h"
#include "mat_elasticity_iso.h"
#include "mesh_mgr.h"
#include "algo_modl.h"
#include "algo_hexvue.h"
#include "algo_refine.h"
#include "algo_show_hier.h"

static char help[] = "Generalized eigenvalue problem Kx=o^2Mx";

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
  int ierr; 

  GCELL_SOLID_H8::register_make_func ();
  ECELL_MODL_H8::register_make_func ();
  //  ECELL_MODL_H8_BBAR::register_make_func ();
  ECELL_HEXVUE_H8<3>::register_make_func ();
  ECELL_HEXVUE_H8<4>::register_make_func ();
  MAT_ELASTICITY_ISO::register_make_func ();

  DB db;
  db.load_file (algorithm + ".fpar");

  MGR mgr (&db);
  SMART_HANDLE<LOGGER_STREAM> toplsb = mgr.logger ("FAMULS (C) 2001-2003 P. Krysl", true);

  ALGO_MODL a (algorithm, &mgr);
  a.setup ();
  a.solve ();
  ALGO_SHOW_HIER sh (a.name(), &mgr);
  sh.set_node_size_mult (node_size_mult); 
  sh.set_fen_rep(FEN::NUMBER); 
  sh.save_as_elixir_file (algorithm + ".egf", a.u()); 

  *toplsb << "done" << endl;
  
  return 0;
}

