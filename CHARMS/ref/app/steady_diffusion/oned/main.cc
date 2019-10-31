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
#include "fles.h"
#include "gmesh.h"
#include "gsubmesh.h"
#include "gcell_line_l2.h"
#include "ecell_steady_diffusion_l2.h"
#include "ecell_errest_steady_diffusion_l2.h"
#include "mgr.h"
#include "algo_steady_diffusion.h"
#include "mat_diffusion_iso.h"
#include "algo_xmgr.h"
#include "proto_field_transfer.h"

static void
xmgr_bfun_set (MGR *mgr, FIELD_SCALAR *phi, FIELD_VECTOR *geometry)
{
  
  FIXED_VECTOR<1>  zero(0.0);
  FIXED_VECTOR<1>  one(1.0);
  for (size_t j = 0; j < phi->npairs(); j++) {
    FIELD_PAIR<1> *pair = phi->field_pair (j);
    pair->set_dofparam (zero);
  }

  ALGO_XMGR xmgr ("", mgr);
  for (size_t j = 0; j < phi->npairs(); j++) {
    FIELD_PAIR<1> *pair = phi->field_pair (j);
    pair->set_dofparam (one);
    char buf[132];
    sprintf (buf, "bfun %d", pair->bfun()->fen()->id());
    xmgr.run (string(buf), phi, geometry);
    pair->set_dofparam (zero);
  }
  
}

#if 0
static void
check_transfers (ALGO_STEADY_DIFFUSION *a, MGR *mgr, ALGO_XMGR *xmgr) 
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
      //cerr << "Transfering " << phi_j->name() << " to " << phi_k->name() << endl;
      //phi_j->debug_display();
      //alt_phi->debug_display();
      t.transfer (phi_j, alt_phi);
      //alt_geometry->debug_display();
      //alt_phi->debug_display();
      xmgr->run (alt_phi->name() + ".xmgr", alt_phi, alt_geometry);
    }
  }
}  
#endif

int
main (int argc, char *argv[])
{
  ARGS::start_process_options (argc, argv, "algorithm name");
  int nadapt = ARGS::get_opt_int ("-nadapt", 3, "# number of adaptation cycles");
  int checkt = ARGS::get_opt_int ("-checkt", 0, "# check transfers");
  ARGS::finish_process_options ();
  
  if (argc < 2) {
    cerr << "Expected algorithm name as the first argument" << endl;
    exit (1);
  }
  string algorithm = argv[1];
  cerr<<"# running"<< algorithm<<endl; 
  fles_initialize (&argc, &argv);

  GCELL_LINE_L2::register_make_func ();
  ECELL_STEADY_DIFFUSION_L2::register_make_func ();
  ECELL_ERREST_STEADY_DIFFUSION_L2::register_make_func ();
  MAT_DIFFUSION_ISO::register_make_func ();

  DB db;
  db.load_file (algorithm + ".fpar");

  MGR mgr (&db);
  SMART_HANDLE<LOGGER_STREAM> toplsb = mgr.logger ("FAMULS (C) 2001-2002 P. Krysl", true);

  ALGO_STEADY_DIFFUSION a (algorithm, &mgr);
 
  a.setup ();
  a.solve ();
  

  ALGO_XMGR xmgr ("", &mgr);
  xmgr.run ("solution", a.phi(), a.geometry());

  string progname = string (argv[0]);
  bool adaptive = false;
  string top = "algorithms/steady_diffusion/" + algorithm;
  string path = top + "/adapt/adaptive";
  if (db.param_defined (path)) {
    adaptive = db.DB_GET_BOOL (path);
  }
  
  if ((progname == "adapt") || (adaptive)) {
    
    //ALGO_REFINE r (a.name(), &mgr, a.gmesh());

    for (int cycle = 0; cycle < nadapt; cycle++) {
      
      SMART_HANDLE<LOGGER_STREAM> lsb = mgr.logger("Adaptation cycle", true);
      *lsb << "cycle " << cycle << endl;
      a.adapt ();
      a.solve ();
      char buf[32]; sprintf (buf, "%d", cycle);
      xmgr.run (string ("adaptive solution ") + string (buf), a.phi(), a.geometry());

      *lsb << "done in adaptation cycle " << cycle << endl;
      
    }
    
  }

  fles_finalize ();
  
  *toplsb << "done" << endl;

  return 0;
}
