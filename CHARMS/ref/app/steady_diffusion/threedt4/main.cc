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
#include "gcell_solid_t4.h"
#include "ecell_steady_diffusion_t4.h"
#include "ecell_errest_steady_diffusion_t4.h"
#include "ecell_hexvue_t4.h"
#include "mgr.h"
#include "algo_steady_diffusion.h"
#include "algo_hexvue.h"
#include "algo_refine.h"
#include "algo_show_hier.h"

typedef CkitPoint3 *PTP;

static int
vdt_4_pt_circumsphere (PTP p0, PTP p1, PTP p2, PTP p3, CkitPoint3 *center)
{
  return G3D_tet_circumsphere (p0, p1, p2, p3, center);
}

static int
vdt_4_pt_insphere (PTP p0, PTP p1, PTP p2, PTP p3, CkitPoint3 *center)
{
  return G3D_tet_insphere (p0, p1, p2, p3, center);
}

static double
vdt_tet_radius_quality (CkitPoint3 pts[4]) 
{
  CkitPoint3 a, b, n= {0,0,0} , v, c;
  double IR, OR;

  if (!vdt_4_pt_circumsphere (&pts[0], &pts[1], &pts[2], &pts[3], &c)) {
    vdt_4_pt_circumsphere (&pts[0], &pts[1], &pts[2], &pts[3], &c);
    cerr << "No circumsphere" << endl;
    return 0;
  } 
  VECMAC_3_VEC (v, =, pts[0], -, c);
  OR = VECMAC_NORM (v);
  if (!vdt_4_pt_insphere (&pts[0], &pts[1], &pts[2], &pts[3], &c)) {
    cerr << "No insphere" << endl;
    return 0;
  } 
  VECMAC_3_VEC (a, =, pts[1], -, pts[0]);
  VECMAC_3_VEC (b, =, pts[2], -, pts[0]);
  VECMAC_CROSS (n, =, a, x, b);
  VECMAC_NORMALIZE (n);
  VECMAC_3_VEC (v, =, c, -, pts[0]);
  IR = VECMAC_DOT (v, n);
  
  return (3 * IR / OR);
} 

static double
vdt_tet_angle_quality_pt (CkitPoint3 *pos1, CkitPoint3 *pos2,
                          CkitPoint3 *pos3, CkitPoint3 *pos4,
                          double *min_angle,
                          double *max_angle,
                          double angles[6])
{
  CkitPoint3 normal={0,0,0}, normal1={0,0,0}, normal2={0,0,0}, normal3={0,0,0};
  CkitPoint3 vec1, vec2, vec3;
  double norm, norm1, norm2, norm3, mina = M_PI, maxa = 0.0;
  double angle1, angle2, angle3, angle4, angle5, angle6;

  *min_angle = 0;
  *max_angle = M_PI;

  /* face normals */
  VECMAC_3_VEC (vec1, =, (*pos1), -, (*pos3));
  VECMAC_3_VEC (vec2, =, (*pos2), -, (*pos3));
  VECMAC_CROSS (normal, =, vec1, x, vec2);
  norm = VECMAC_NORM (normal);
  if (norm <= 0) {
    cerr << "Zero-area tet face" << endl;
    return 0;
  }
  VECMAC_3_VEC (vec1, =, (*pos1), -, (*pos4));
  VECMAC_3_VEC (vec2, =, (*pos2), -, (*pos4));
  VECMAC_3_VEC (vec3, =, (*pos3), -, (*pos4));
  VECMAC_CROSS (normal1, =, vec3, x, vec2);
  norm1 = VECMAC_NORM (normal1);
  VECMAC_CROSS (normal2, =, vec1, x, vec3);
  norm2 = VECMAC_NORM (normal2);
  VECMAC_CROSS (normal3, =, vec2, x, vec1);
  norm3 = VECMAC_NORM (normal3);
  if (norm1 <= 0 || norm2 <= 0 || norm3 <= 0) {
    cerr << "Zero-area tet face" << endl;
    return 0;
  }

  /* convert to unit normals */
  VECMAC_SCALAR (normal, /=, norm);
  VECMAC_SCALAR (normal1, /=, norm1);
  VECMAC_SCALAR (normal2, /=, norm2);
  VECMAC_SCALAR (normal3, /=, norm3);

  /* compute angles */
  angle1 = VECMAC_DOT (normal, normal1);
  angle2 = VECMAC_DOT (normal, normal2);
  angle3 = VECMAC_DOT (normal, normal3);
  angle4 = VECMAC_DOT (normal2, normal3);
  angle5 = VECMAC_DOT (normal3, normal1);
  angle6 = VECMAC_DOT (normal1, normal2);

#define EVA(N)                                                          \
  {                                                                     \
    if (fabs (angle##N) > 1) angle##N = (angle##N > 0 ? +1 : -1);       \
    angle##N = acos (-angle##N);                                        \
  }
  EVA (1); EVA (2); EVA (3); EVA (4); EVA (5); EVA (6);

  /* evaluate min/max */
  mina = min (mina, angle1);
  mina = min (mina, angle2);
  mina = min (mina, angle3);
  mina = min (mina, angle4);
  mina = min (mina, angle5);
  mina = min (mina, angle6);

  maxa = max (maxa, angle1);
  maxa = max (maxa, angle2);
  maxa = max (maxa, angle3);
  maxa = max (maxa, angle4);
  maxa = max (maxa, angle5);
  maxa = max (maxa, angle6);

  /* output arguments */
  *min_angle = mina;
  *max_angle = maxa;
  angles[0] = angle1;
  angles[1] = angle2;
  angles[2] = angle3;
  angles[3] = angle4;
  angles[4] = angle5;
  angles[5] = angle6;

  return (maxa > 0 ? mina / maxa : 0);
}

static void
mesh_quality_report (GMESH *gmesh)
{
  int n = 0;
  GMESH::gsubmesh_enumerator_t si = gmesh->gsubmesh_enumerator ();
  GSUBMESH *gsubmesh;
  si.reset();
  while (gsubmesh = si.next()) {
    GSUBMESH::gcell_group_enumerator_t ggi = gsubmesh->gcell_group_enumerator ();
    GCELL_GROUP *gcell_group;
    ggi.reset();
    while (gcell_group = ggi.next()) {
      GCELL_GROUP::gcell_enumerator_t gi = gcell_group->gcell_enumerator ();
      GCELL_SOLID_T4 *gcell;
      gi.reset();
      while (gcell = dynamic_cast <GCELL_SOLID_T4 *> (gi.next())) {
        const CONN_BASE *conn = gcell->conn();
        CkitPoint3 p[4];
        for (int i = 0; i < 4; i++) {
          FEN *fen = conn->fen(i);
          POINT ref_loc = fen->ref_loc();
          p[i].x = ref_loc(0);
          p[i].y = ref_loc(1);
          p[i].z = ref_loc(2);
        }
        double radqual = vdt_tet_radius_quality (p);
        cout << n++ << " " << radqual << endl;
      }
    }
  }
}

int
main (int argc, char **argv)
{
  ARGS::start_process_options (argc, argv, "algorithm name");
  int nadapt = ARGS::get_opt_int ("-nadapt", 3, "# number of adaptation cycles");
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

  GCELL_SOLID_T4::register_make_func ();
  ECELL_STEADY_DIFFUSION_T4::register_make_func ();
  ECELL_HEXVUE_T4<1>::register_make_func ();
  ECELL_ERREST_STEADY_DIFFUSION_T4::register_make_func ();
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
    
    // ALGO_REFINE r (a.name(), &mgr, a.gmesh());
    ALGO_SHOW_HIER sh (a.name(), &mgr);

    sh.set_node_size_mult (node_size_mult);
    sh.set_level_as_fen_layer (level_as_fen_layer);
    sh.set_fen_rep (static_cast<FEN::FEN_REP> (fen_rep));
    sh.save_as_elixir_file (algorithm + "-hier.egf", a.phi());
    
    for (int cycle = 0; cycle < nadapt; cycle++) {
      char buf[64];
      
      SMART_HANDLE<LOGGER_STREAM> lsb = mgr.logger("Adaptation cycle", true);
      *lsb << "cycle " << cycle << endl;
      
      sprintf (buf, "-meshsize%d.hvu", cycle);
      a.adapt ();

      sprintf (buf, "-hier%d.egf", cycle);
      sh.set_node_size_mult (node_size_mult);
      sh.set_level_as_fen_layer (level_as_fen_layer);
      sh.set_fen_rep (static_cast<FEN::FEN_REP> (fen_rep));
      sh.save_as_elixir_file (algorithm + string (buf), a.phi());

      a.solve ();
      
      sprintf (buf, "-adapt%d.hvu", cycle);
      h.write_field_to_name (algorithm + string (buf), a.geometry(), a.geometry(), a.phi());
      
      *lsb << "done in adaptation cycle " << cycle << endl;
      
    }
    
  }

  mesh_quality_report(a.gmesh());

  fles_finalize ();

  *toplsb << "done" << endl;

  return 0;
}
