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
#include "gcell_surf_q4.h"
#include "ecell_steady_diffusion_q4.h"
#include "ecell_errest_steady_diffusion_q4.h"
#include "ecell_hexvue_q4.h"
#include "mgr.h"
#include "algo_steady_diffusion.h"
#include "algo_hexvue.h"
#include "algo_refine.h"
#include "algo_show_hier.h"
#include "proto_field_transfer.h"

#if 0
static void
check_transfers (ALGO_STEADY_DIFFUSION *a, MGR *mgr) 
{
  cerr << "Checking transfers for " << a->nphi() << " fields" << endl;
  for (size_t j = 0; j < a->nphi(); j++) {
    FIELD_SCALAR *phi_j = a->phi(j);
    for (size_t k = 0; k < a->nphi(); k++) {
      char buf[512];
      sprintf (buf, "phi%d_to_phi%d", j, k);
      PROTO_FIELD_TRANSFER<1> t;
      FIELD_SCALAR *phi_k = a->phi(k);
      FIELD_SCALAR *alt_phi = phi_k->clone(string(buf));
      alt_phi->setto (0);
      FIELD_VECTOR *alt_geometry = a->geometry(k);
      t.transfer (phi_j, alt_phi);
      ALGO_HEXVUE<1> h ("hexvue", mgr);
      alt_geometry->debug_display();
      alt_phi->debug_display();
      h.write_field_to_name (alt_phi->name() + ".hvu", alt_geometry, alt_geometry, alt_phi);
    }
  }
}  
#endif

int
main (int argc, char **argv)
{
  ARGS::start_process_options (argc, argv, "fpar_file");
  int nadapt = ARGS::get_opt_int ("-nadapt", 3, "# number of adaptation cycles");
  ARGS::finish_process_options ();
  
  if (argc < 2) {
    cerr << "Expected algorithm name as the first argument" << endl;
    exit (1);
  }
  string algorithm = argv[1];

  fles_initialize (&argc, &argv);

  GCELL_SURF_Q4::register_make_func ();
  ECELL_STEADY_DIFFUSION_Q4::register_make_func ();
  ECELL_HEXVUE_Q4<1>::register_make_func ();
  ECELL_ERREST_STEADY_DIFFUSION_Q4::register_make_func ();
  MAT_DIFFUSION_ISO::register_make_func ();

  DB db;
  db.load_file (algorithm + ".fpar");

  MGR mgr (&db);
  SMART_HANDLE<LOGGER_STREAM> toplsb = mgr.logger ("FAMULS (C) 2001-2002 P. Krysl", true);

  ALGO_STEADY_DIFFUSION a (algorithm, &mgr);
  a.setup ();
  a.solve ();
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
    
    // ALGO_REFINE r (a.name(), &mgr, a.gmesh()); //RAW
    ALGO_SHOW_HIER sh (a.name(), &mgr);
    sh.set_fen_rep(FEN::NUMBER); 
    for (int cycle = 0; cycle < nadapt; cycle++) {
      char buf[32];
      SMART_HANDLE<LOGGER_STREAM> lsb = mgr.logger("Adaptation cycle", true);
      *lsb << "cycle " << cycle << endl;
      a.adapt (); //(&r) :RAW
      a.solve ();
      sprintf (buf, "-adapt%d.hvu", cycle);
      h.write_field_to_name (algorithm + string (buf), a.geometry(), a.geometry(), a.phi());
      sprintf (buf, "-hier%d.egf", cycle);
      sh.save_as_elixir_file (algorithm + string (buf), a.phi());

      *lsb << "done in adaptation cycle " << cycle << endl;
      
    }
    
  }

  fles_finalize ();

  *toplsb << "done" << endl;

  return 0;
}
