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
#include "mgr.h"
#include "algo_steady_diffusion.h"
#include "gcell_surf_t6.h"
#include "ecell_steady_diffusion_t6.h"
#include "ecell_hexvue_t6.h"
#include "ecell_errest_steady_diffusion_t6.h"
#include "algo_hexvue.h"
#include "algo_refine.h"
#include "algo_show_hier.h"

int
main (int argc, char **argv)
{
  ARGS::start_process_options (argc, argv, "fpar_file");
  int nadapt = ARGS::get_opt_int ("-nadapt", 3, "# number of adaptation cycles");
  int whier = ARGS::get_opt_int ("-whier", 0, "# write mesh hierarchy; true or false");
  bool level_as_fen_layer = ARGS::get_opt_bool ("-level_as_fen_layer", true, "should fen be placed in layers?");
  double node_size_mult = ARGS::get_opt_double ("-node_size_mult", 0.06, "node size multiplier");
  int fen_rep = ARGS::get_opt_int ("-fen_rep", 0, "FEN representation (0=marker,1=ball,2=number)");
  ARGS::finish_process_options ();
  
  if (argc < 2) {
    cerr << "Expected algorithm name as the first argument" << endl;
    exit (1);
  }
  string algorithm = argv[1];

  fles_initialize (&argc, &argv);

  GCELL_SURF_T6::register_make_func ();
  ECELL_STEADY_DIFFUSION_T6::register_make_func ();
  ECELL_HEXVUE_T6<1>::register_make_func ();
  ECELL_ERREST_STEADY_DIFFUSION_T6::register_make_func ();
  MAT_DIFFUSION_ISO::register_make_func ();

  DB db;
  db.load_file (algorithm + ".fpar");

  MGR mgr (&db);
  SMART_HANDLE<LOGGER_STREAM> toplsb = mgr.logger ("FAMULS (C) 2001-2002 P. Krysl", true);

  ALGO_STEADY_DIFFUSION a (algorithm, &mgr);
  a.setup ();
  a.solve ();
  {
    ALGO_SHOW_HIER sh (a.name(), &mgr);
    sh.set_node_size_mult (node_size_mult);
    sh.set_level_as_fen_layer (level_as_fen_layer);
    sh.set_fen_rep (static_cast<FEN::FEN_REP> (fen_rep));
    sh.save_as_elixir_file (algorithm + ".egf", a.phi());
  }
  ALGO_HEXVUE<1> h ("hexvue", &mgr);
  h.write_field (a.geometry (), a.geometry (), a.phi ());

  string progname = string (argv[0]);
  bool adaptive = false;
  string top = "algorithms/steady_diffusion/" + algorithm;
  string path = top + "/adapt/adaptive";
  if (db.param_defined (path)) {
    adaptive = db.DB_GET_BOOL (path);
  }
  
  if ((progname == "adapt") || (adaptive)) {
    
  //  ALGO_REFINE r (a.name(), &mgr, a.gmesh());
    ALGO_SHOW_HIER sh (a.name(), &mgr);

    for (int cycle = 0; cycle < nadapt; cycle++) {
      char buf[32];
      
      SMART_HANDLE<LOGGER_STREAM> lsb = mgr.logger("Adaptation cycle", true);
      *lsb << "cycle " << cycle << endl;
      
      a.adapt ();
      a.solve ();
      
      sprintf (buf, "-adapt%d.hvu", cycle);
      h.write_field_to_name (algorithm + string (buf), a.geometry(), a.geometry(), a.phi());
      
      sprintf (buf, "-hier%d.egf", cycle);
      sh.set_node_size_mult (node_size_mult);
      sh.set_level_as_fen_layer (level_as_fen_layer);
      sh.set_fen_rep (static_cast<FEN::FEN_REP> (fen_rep));
      sh.save_as_elixir_file (algorithm + string (buf), a.phi());

      *lsb << "done in adaptation cycle " << cycle << endl;
      
    }
    
  }

  fles_finalize ();

  *toplsb << "done" << endl;

  return 0;
}
