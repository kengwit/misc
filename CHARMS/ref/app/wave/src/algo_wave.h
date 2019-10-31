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
#ifndef ALGO_WAVE_H
# define ALGO_WAVE_H

#include <list>
#include <iterator>
#include "field.h"
#include "gcell.h"
#include "proto_wave.h"
#include "algo.h"
#include "db.h"
#include "ebc.h"
#include "algo_refine.h"
#include "mesh_mgr.h"
#include "watchpoint_ecell.h"
#include "algo_xhexvue.h"
#include "algo_motion_pict.h"

class ALGO_WAVE : public ALGO {

 public: // object functions ////////////////////////////////////////

  /**
     Construct an algorithm.  Give it a name, and let it load
     its parameters from the database.  The mesh given on input
     must be the one specified for the algorithm in the database: an
     exception is raised if that is not the case.
  */
  ALGO_WAVE (string name, MGR *mgr);
  /**
   */
  ~ALGO_WAVE ();
  /**
     Prepare the algorithm for solution.  
  */
  void setup ();
  /**
     Solve for the primary field.
  */
  void solve (double tstart, double tend);

  /**
   solve and adapt the solution
  */
  void run(); 
  /**
     Return the geometry.
  */
  FIELD_VECTOR * geometry () { return _geometry; }
  /**
     Return the primary field (displacements).
  */
  FIELD_VECTOR *u () { return _u; }
  /**
     Return the mesh.
  */
  GMESH *gmesh () { return _gmesh; }
  /**
     Adapt the basis function sets for the primary field and the geometry.
  */
  void adapt ();
  /**
     Apply initial condition on displacement and velocity field 
   */
  void apply_initial_cond ();
  /**
    Get the time step.
   */
  double dt ();
  /** Get current time.
   */
  double t() { return _t;}
  /**
     Protocol may use this.
  */
  void  output();
  /**
     Protocol may use this.
  */
  void  increment_t();
  /**
   */
  SMART_HANDLE<EBC<3> > ebc() { return _ebc; }
  /**
    fraction that specifies the time step as a percentage of stable time.  
   */
  double red_fact(); 
  double tstart() {return _tstart;}
  sizet buffer_size () {return _buffer_size;}
  void   set_nadapt(sizet nadapt);
  void set_hist_file_prefix (string hist_file_prefix);

 private: // object functions /////////////////////////////////////

  vector<double>  make_tlist();
  void update_history();
  void flush_history();
  void map_watchpoints ();
  void unmap_watchpoints();
  void update_geometry();
  void assemble_error();

 private: // object data //////////////////////////////////////////

  GMESH                     *_gmesh;
  list <GSUBMESH *>          _gsubmeshes;
  FIELD_VECTOR              *_geometry;
  FIELD_VECTOR              *_u;
  FIELD_VECTOR              *_v;
  SMART_HANDLE<EBC<3> >      _ebc;
  ALGO_ERREST<3>            *_algo_errest;
  ALGO_REFINE               *_algo_refine;
  PROTO_WAVE                *_proto_wave;
  double                     _dt; 
  double                     _maxdt; 
  double                     _tstart; 
  double                     _t;
  double                     _tend;   
  double                     _red_fact; 
  FUNC<3>                   *_u0;
  FUNC<3>                   *_v0;
  list<WATCH_POINT_ECELL*>   _watch_point;
  WATCH_POINT               *_global_watch_point; 
  string                     _hist_file_prefix;
  double                     _time_between_outputs;
  double                     _time_between_graphics;
  int                        _n_output_times;
  int                        _n_graphics_times;
  sizet                      _buffer_size; 
  double                     _prescribed_dt;
  sizet                      _nadapt;
  sizet                      _ncycle;
  double                     _back_track_frac;
  vector<double>             _tlist;
  bool                       _is_back_tracking;
  bool                       _force_back_track_output;
  ALGO_MOTION_PICT          *_algo_motion_pict;
  bool                       _motion_pict_show;
  bool                       _motion_pict_stop_in_solve;
  FIELD_VECTOR              *_updated_geometry;
  bool                       _accumulate_error;

};

#endif
