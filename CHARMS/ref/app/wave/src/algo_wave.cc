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
#include "algo_wave.h"
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
#include "tokensP.h"
#include "algo_hexvue.h"
#include "ecell_hexvue_h8.h"
#include "algo_show_hier.h"

ALGO_WAVE::ALGO_WAVE (string name, MGR *mgr) : ALGO (name, mgr)
{
  CHECK (this->mgr()->db() != 0, EXCEPTION_NULL_PTR,;);
  _gmesh = 0;
  string top = "algorithms/wave/" + this->name ();
  string gmesh_name = this->mgr()->db()->DB_GET_STRING (top + "/gmesh");
  _gmesh = this->mgr()->mesh_mgr()->gmesh (gmesh_name);
  _gsubmeshes.clear();
  _algo_refine = new ALGO_REFINE(name, mgr, _gmesh);
  _geometry = 0;
  _u = 0;
  _v = 0;
  _ebc = 0;
  _algo_errest = 0;
  _t           = 0;
  _dt          = 0;
  _tstart      = 0;
  _tend        = 0;
  _ncycle      = 0;
  _red_fact    = 1;
  _u0 = 0;
  _v0 = 0;
  _watch_point.clear();
  _hist_file_prefix = "";
  _time_between_outputs = 0;
  _time_between_graphics = 0;
  _buffer_size = 100;
  _nadapt = 0;
  _ncycle = 0;
  _is_back_tracking = false;
  _tlist.clear();
  _algo_motion_pict = 0;
  _motion_pict_show = false;
  _motion_pict_stop_in_solve = false;
  _updated_geometry = 0;
  _accumulate_error = false;
}

void
ALGO_WAVE::setup ()
{
  CHECK (this->mgr()->db() != 0, EXCEPTION_NULL_PTR,;);
  DB *db = this->mgr()->db();
  string top = "algorithms/wave/" + name ();
  
  // construct the gsubmesh list
  list <string> gsubmesh_names = db->DB_GET_STRING_LIST (top + "/gsubmeshes");
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
  SMART_HANDLE<BFUN_SET> geom_bfun_set
    = BFUN_SET::make (_algo_refine->refinement_strategy(),_gmesh);
  _geometry = new FIELD_VECTOR ("geometry", geom_bfun_set, zero);
  ALGO_GEOM init_geom ("", this->mgr());
  init_geom.init (_geometry, _gmesh);
  
  // initial condition functions
  string expr[3];  
  list<string> u0_  = db->DB_GET_STRING_LIST (top + "/initial_conditions/u0");
  list<string>::iterator i = u0_.begin();  
  expr[0] = (*i);i++; 
  expr[1] = (*i);i++;
  expr[2] = (*i);
  _u0 = new FUNC<3> (expr);
  
  list<string> v0_  = db->DB_GET_STRING_LIST (top + "/initial_conditions/v0");
  list<string>::iterator i1 = v0_.begin();  
  expr[0] = (*i1);i1++; 
  expr[1] = (*i1);i1++;
  expr[2] = (*i1);
  _v0 = new FUNC<3> (expr);
  
  // create the fields u and v
  SMART_HANDLE<BFUN_SET> u_bfun_set
    = BFUN_SET::make (_algo_refine->refinement_strategy(),_gsubmeshes);
  _u = new FIELD_VECTOR ("u", u_bfun_set, zero);
  _v = new FIELD_VECTOR ("v", u_bfun_set, zero);
  
  // read watchpoint data
  _watch_point.clear();
  if (db->param_defined (top+ "/watchpoints/list")) {  
    list<string> watchpoints =  db->DB_GET_STRING_LIST (top+ "/watchpoints/list");
    _time_between_outputs = 0;
    _time_between_outputs = db->DB_GET_DOUBLE (top+ "/watchpoints/time_between_outputs");
    _time_between_graphics = 0;
    _time_between_graphics = db->DB_GET_DOUBLE (top+ "/watchpoints/time_between_graphics");
    _buffer_size = 100;
    _buffer_size = db->DB_GET_INTEGER (top+ "/watchpoints/buffer_size");
    _force_back_track_output = false;
    _force_back_track_output = db->DB_GET_BOOL (top+ "/watchpoints/force_back_track_output");
    if (db->param_defined (top + "/watchpoints/hist_file_prefix")) {
      _hist_file_prefix = db->DB_GET_STRING (top+ "/watchpoints/hist_file_prefix");
      _hist_file_prefix = _hist_file_prefix + "_";
    }
    list<string> global_watch_point_list;
    global_watch_point_list.push_back("kinetic_energy"); global_watch_point_list.push_back("strain_energy"); 
    global_watch_point_list.push_back("total_energy");   global_watch_point_list.push_back("momentum_x");
    global_watch_point_list.push_back("momentum_y");  global_watch_point_list.push_back("momentum_z");
    global_watch_point_list.push_back("nbfuns");
    global_watch_point_list.push_back("mass");
    _global_watch_point = new WATCH_POINT  (global_watch_point_list, _hist_file_prefix+ name() );
    tokens_parser parser = tokens_new_parser ();
    vector <sizet> gcell_ids;
    for (list<string>::iterator i = watchpoints.begin(); i != watchpoints.end(); i++) {
      tokens_parse_line (parser, (char *) (*i).c_str ());
      CHECK (tokens_total_of_tokens (parser) == 1, EXCEPTION_BAD_VALUE,;);
      string wp =   tokens_token_as_string(parser,1);
      list <string>var_list =   db->DB_GET_STRING_LIST (top+"/watchpoints/"+wp+"/var");
      list<string> correct_var_list;
      for (list<string>::iterator ix = var_list.begin(); ix != var_list.end(); ix++) {
        tokens_parse_line (parser, (char *) (*ix).c_str ());
        string w  =   tokens_token_as_string(parser,1);
        correct_var_list.push_back(w);
      }
      sizet gcell_id =0;
      if (db->param_defined(top+"/watchpoints/"+wp+"/gcell_id")) {
        gcell_id = (sizet) db->DB_GET_DOUBLE (top+"/watchpoints/"+wp+"/gcell_id");
        gcell_ids.push_back(gcell_id);
      }
      string prefix = _hist_file_prefix + name() + "_" + wp;
      if (db->param_defined(top+"/watchpoints/"+wp+"/param_loc")) {
        list <double>param_locs = db->DB_GET_DOUBLE_LIST (top+"/watchpoints/"+wp+"/param_loc");
        POINT param_loc;
        int j = 0; 
        for (list<double>::iterator i1 = param_locs.begin(); i1!= param_locs.end(); i1++) {
          param_loc(j) = (*i1);
          j ++;
        }
        WATCH_POINT_ECELL *watch_point
          = new WATCH_POINT_ECELL (correct_var_list, prefix, gcell_id, param_loc);
        watch_point->set_buf_size(_buffer_size);
        _watch_point.push_back(watch_point);
      } else if (db->param_defined(top+"/watchpoints/"+wp+"/coord")){
        list <double>coords = db->DB_GET_DOUBLE_LIST (top+"/watchpoints/"+wp+"/coord");
        list <GSUBMESH *>::const_iterator si = _gsubmeshes.begin ();
        int counter = 0;
        GCELL *fen_gcell = 0;
        FEN *mapped_fen = 0;
        while (si != _gsubmeshes.end ()) {
          GSUBMESH *gsubmesh = *si;
          GSUBMESH::gcell_group_enumerator_t e = gsubmesh->gcell_group_enumerator ();
          GCELL_GROUP *gcell_group;
          e.reset ();
          while ((gcell_group = e.next ())) {
            GCELL_GROUP::gcell_enumerator_t e1 = gcell_group->gcell_enumerator ();
            GCELL *gcell;
            e1.reset ();
            while ((gcell = e1.next ())) {
              counter++;
              for (sizet gcell_conn = 0; gcell_conn != gcell->conn()->nfens(); gcell_conn++) {
                FEN *fenx = gcell->conn()->fen(gcell_conn); 
                POINT loc = fenx->ref_loc();
                bool find_fen = true;
                int temp_int = 0;
                for (list<double>::iterator i2 = coords.begin(); i2 != coords.end(); i2++) {
                  if ( loc(temp_int) != *i2) find_fen = false; temp_int ++;    
                }
                if (find_fen) {
                  fen_gcell = gcell;
                  mapped_fen = fenx;
                  goto exit;
                }
              }
            }
          }
        }
      exit:
        POINT param_loc;
        fen_gcell->map_fen(mapped_fen, &param_loc);
        sizet gcell_id = (sizet)fen_gcell;
      back: 
        for (vector<sizet>::iterator gi = gcell_ids.begin(); gi != gcell_ids.end(); gi++) {
          if (gcell_id == (*gi)) {
            gcell_id = gcell_id+1;
            goto back;
          }
        }
        gcell_ids.push_back(gcell_id);
        fen_gcell->set_id(gcell_id);
        WATCH_POINT_ECELL *watch_point
          = new WATCH_POINT_ECELL (correct_var_list, prefix, gcell_id, param_loc);
        watch_point->set_buf_size(_buffer_size);
        _watch_point.push_back(watch_point);    
      } else {
        cout<<"Error:Neither parametric location  nor position coordinate defined for watchpoint"<<endl;
        break;
      }
    }
    tokens_delete_parser(parser);
  }
  
  //apply initial conditions
  apply_initial_cond();
   
  // read the essential boundary conditions
  string path = top + "/essential_bcs/ebc_file";
  string file = db->DB_GET_STRING (path);
  _ebc = new EBC<3> (_gmesh);
  CHECK (_ebc->read (file), EXCEPTION_BAD_VALUE,;);

  _nadapt = 0;
  if (db->param_defined (top + "/time_step/nadapt")) {  
    _nadapt = db->DB_GET_INTEGER (top + "/time_step/nadapt");
  }
 
  _prescribed_dt = 0;
  if (db->param_defined (top + "/time_step/dt")) {  
    _prescribed_dt = db->DB_GET_DOUBLE (top + "/time_step/dt");
  }
 
  _maxdt = 0;
  if (db->param_defined (top + "/time_step/maxdt")) {  
    _maxdt = db->DB_GET_DOUBLE (top + "/time_step/maxdt");
  }
 
  _tstart = 0;
  if (db->param_defined (top + "/time_step/tstart")) { 
    _tstart = db->DB_GET_DOUBLE (top + "/time_step/tstart");
  }

  _tend = 1;
  if (db->param_defined (top + "/time_step/tend")) { 
    _tend = db->DB_GET_DOUBLE (top + "/time_step/tend");
  }

  _red_fact = 1;
  if (db->param_defined (top + "/time_step/red_fact")) { 
    _red_fact = db->DB_GET_DOUBLE (top + "/time_step/red_fact");
  }
   
  _back_track_frac = 0;
  if (db->param_defined (top + "/time_step/back_track_frac")) { 
    _back_track_frac = db->DB_GET_DOUBLE (top + "/time_step/back_track_frac");
  }

  _tlist.clear();
  if (db->param_defined (top + "/time_step/tlist")) { 
    list<double> tl= db->DB_GET_DOUBLE_LIST (top + "/time_step/tlist");
    {    list<double>::iterator i = tl.begin();
      if (*i != _tstart) _tlist.push_back(_tstart);// make sure the list starts w/ _tstart
    }
    double pt=-FLT_MAX;
    for (list<double>::iterator i = tl.begin(); i != tl.end();  i++) {
      CHECK((*i) > pt, EXCEPTION_BAD_VALUE, );
      _tlist.push_back(*i);
      pt = *i;
    }
    {    list<double>::iterator i = tl.end();
      if (*i != _tend) _tlist.push_back(_tend);// make sure the list ends w/ _tend
    }
    _nadapt = _tlist.size();
  }

  _motion_pict_show = false;
  if (db->param_defined (top + "/motion_pict/show")) { 
    _motion_pict_show = db->DB_GET_BOOL (top + "/motion_pict/show");
  }
  _motion_pict_stop_in_solve = false;
  if (db->param_defined (top + "/motion_pict/stop_in_solve")) { 
    _motion_pict_stop_in_solve = db->DB_GET_BOOL (top + "/motion_pict/stop_in_solve");
  }

  if (_motion_pict_show) {
    _algo_motion_pict = new ALGO_MOTION_PICT (this->name(), this->mgr(), this->_gmesh);
    if (db->param_defined (top + "/motion_pict/scale")) { 
      double scale = db->DB_GET_DOUBLE (top + "/motion_pict/scale");
      _algo_motion_pict->set_scale(scale);
    }
  }

  {
    string path = "algorithms/refine/" + name();
    if (this->mgr()->db()->param_defined (path + "/accumulate_error")) { 
      _accumulate_error = this->mgr()->db()->DB_GET_BOOL (path + "/accumulate_error");
    }
  }


  _n_graphics_times = 0; 
  _n_output_times = 0;                             

  _t = _tstart;
  _proto_wave = 0;

  if (_motion_pict_show) {
    _algo_motion_pict->event_loop(true, "Ready to run in ALGO_WAVE::setup(); ctrl-p to proceed");
  }
}

void
ALGO_WAVE::update_geometry() {
  if (_updated_geometry == 0) _updated_geometry = _geometry->clone("updated_geometry");
  for (sizet j=0; j<_u->npairs();j++) {
    FIXED_VECTOR<3> ut = _u->field_pair(j)->dofparam();
    FIXED_VECTOR<3> x = _geometry->field_pair(j)->dofparam();
    x.add(_algo_motion_pict->scale(), ut);
    _updated_geometry->field_pair(j)->set_dofparam(x);
  }
}

void
ALGO_WAVE::solve (double tstart, double tend)
{
  CHECK (this->mgr()->db() != 0, EXCEPTION_NULL_PTR,;);
  string top = "algorithms/wave/" + name ();

  SMART_HANDLE<LOGGER_STREAM> lsb = this->mgr()->logger("solve", true);
  
  if (_proto_wave) delete _proto_wave;
  
  // Set up protocol  
   _proto_wave = new PROTO_WAVE (this, this->mgr()->db(), _geometry, _u, _v);

  // Now time step
  *lsb << "start=" << tstart << "; end=" << tend << endl_notime;
  if ( _prescribed_dt !=0)_dt = _prescribed_dt;
  else _dt = (_proto_wave->stable_time_step())*(_red_fact);
  if ((_maxdt > 0 ) && (_dt > _maxdt)) _dt = std::min(_dt,_maxdt); // clamp if needed
  sizet nsteps = (sizet) ((tend - tstart)/this->dt()); 
  *lsb << "time step= " << _dt << " (" << nsteps << " steps)" << endl_notime;
  
  // Apply EBC to the unknown field
  *lsb << "applying EBC's ... ";
  _ebc->apply (_u, _t); 
  _ebc->apply (_v, _t); // RAW have to treat velocities consistently
   
  //_u->attach_ebc(_ebc);
  //_v->attach_ebc(_ebc);

  *lsb << "done" << endl;
  
   if (_algo_motion_pict) {
     update_geometry ();
     _proto_wave->add_graphics (_algo_motion_pict, _updated_geometry);
     char buf[1024];
     sprintf(buf, "Ready to run in ALGO_WAVE::solve() from %g to %g; ctrl-p to proceed", tstart, tend);
     _algo_motion_pict->event_loop((tstart == _tstart) || _motion_pict_stop_in_solve, buf);
   }

  // map  watch_points to ecells
  map_watchpoints();

  // advance the solution
  CHECK(_proto_wave->advance (nsteps), EXCEPTION_BAD_VALUE,;);

  unmap_watchpoints ();
 
  *lsb << "done"  << endl;
}

void 
ALGO_WAVE::map_watchpoints ()
{
  _proto_wave->map_watchpoints(_watch_point);
}

void 
ALGO_WAVE::unmap_watchpoints ()
{
  for (list<WATCH_POINT_ECELL*>::iterator i= _watch_point.begin(); i != _watch_point.end(); i++) {
    (*i)->set_ecell(0);
  }  
}

void
ALGO_WAVE::assemble_error ()
{
  if (_algo_errest == 0) {
    string path = "algorithms/refine/" + this->name();
    string errestim = "gradation";
    if (this->mgr()->db()->param_defined (path + "/errestim")) {
      errestim = this->mgr()->db()->DB_GET_STRING (path + "/errestim");
    }
    if (errestim == "gradation") {
      _algo_errest = new ALGO_ERREST_GRADATION<3> (this->name(), this->mgr(), _gmesh);
    } else {
      _algo_errest = new ALGO_ERREST_ENERGY<3> (this->name(), this->mgr());
    }
  }
  
  // Assemble the error, and then use it to update the refinement context
  if (! _accumulate_error) _algo_errest->assemble_error (_geometry, _u, _t);
  else                     _algo_errest->accumulate_error (_geometry, _u, _t);
}

 
void
ALGO_WAVE::adapt ()
{
  SMART_HANDLE<LOGGER_STREAM> lsb = this->mgr()->logger("adapt", true);
  if (! _accumulate_error) {
    *lsb << "assembling error ... ";
    assemble_error();
    *lsb << "done" << endl;
  }
  *lsb << "updating refinement context ... ";

#if 0
  ALGO_HEXVUE<1> h ("hexvue", this->mgr());
  SMART_HANDLE<FIELD<1> > hbar = _algo_errest->hbar(_u->bfun_set(), 4000);
  h.write_field_to_name ("hexvue", _geometry, _geometry , hbar.get_ptr());
#endif  

  SMART_HANDLE <BFUN_SET> new_bfun_set = _algo_refine->adapt(_algo_errest, _u); 
  *lsb << "done" << endl;
  delete _algo_errest; _algo_errest = 0;

  *lsb << "making new fields ... ";
  FIXED_VECTOR<3>  zero(0);
  FIELD_VECTOR *new_u = new FIELD_VECTOR ("u", new_bfun_set, zero);
  FIELD_VECTOR *new_v = new FIELD_VECTOR ("v", new_bfun_set, zero);
  FIELD_VECTOR *new_geometry = new FIELD_VECTOR ("geom", new_bfun_set, zero);
  *lsb << "done" << endl;
    
  PROTO_FIELD_TRANSFER<3> transf3;
  *lsb << "transfering geometry ... ";
  transf3.transfer (_geometry, new_geometry);
  *lsb << "done" << endl;
  *lsb << "transfering u ... ";
  transf3.transfer (_u, new_u);
  *lsb << "done" << endl;
  *lsb << "transfering v ... ";
  transf3.transfer (_v, new_v);
  *lsb << "done" << endl;
  
   if (_algo_motion_pict) {
     _algo_motion_pict->clear();
   }

  delete _u;
  _u = new_u;
  delete _v;
  _v = new_v;
  delete _geometry;
  _geometry = new_geometry;
  if (_updated_geometry) {
    delete _updated_geometry;
  }
  _updated_geometry = 0;
  
  *lsb << "done" << endl;
}



void ALGO_WAVE::apply_initial_cond ()
{
  for (sizet j=0;j<_u->npairs();j++) {
    FIELD_PAIR<3> *fp = _u->field_pair(j);
    POINT at = fp->bfun()->fen()->ref_loc();
    fp->set_dofparam ((*_u0)(at));
  }
  for (sizet j=0;j<_v->npairs();j++) {
    FIELD_PAIR<3> *fp = _v->field_pair(j);
    POINT at = fp->bfun()->fen()->ref_loc();
    fp->set_dofparam ((*_v0)(at));
  }  
}


double  ALGO_WAVE::dt () {
 return _dt;
} 

double  ALGO_WAVE::red_fact () {return _red_fact;} 


vector<double> ALGO_WAVE::make_tlist()
{
  if (_tlist.empty()) {
    vector<double> tl;
    tl.push_back(_tstart);
    for (sizet i = 1; i <= _nadapt; i++) {
      tl.push_back(_tstart + ((_tend-_tstart)*i/_nadapt));
    }
    return tl;
  } else {
    _nadapt = _tlist.size();
    return _tlist;
  }
}

void ALGO_WAVE::run() 
{
  vector<double> tl = make_tlist();
  SMART_HANDLE<LOGGER_STREAM> lsb = this->mgr()->logger("run", true);
  *lsb << _nadapt <<  " adaptation steps" << endl_notime;
  if (_nadapt==0) solve(_tstart, _tend);
  else {
    if (_back_track_frac != 0) {
      for (sizet i = 1; i <= _nadapt; i++) {
        double tb = tl[i-1];
        double te = tl[i];
        double ti = te - tb;  
        double btt = _back_track_frac*ti;
        FIELD_VECTOR *restart_u = _u->clone("u");
        FIELD_VECTOR *restart_v = _v->clone("v");
        _is_back_tracking = true; //true for no overwriting   
        *lsb << "with current mesh" << endl_notime;
        _t = tb;   
        solve (tb,tb+btt);  
        adapt(); 
        PROTO_FIELD_TRANSFER<3> transf3; 
        transf3.transfer (restart_u, _u); 
        transf3.transfer (restart_v, _v); 
        delete restart_u; delete restart_v;
        _is_back_tracking = false; 
        *lsb << "with adapted mesh" << endl_notime;
        _t = tb;   
        solve(tb,te);
        tb = te;
      }   
    } else {
      for (sizet i = 1; i <= _nadapt; i++) {
        double tb = tl[i-1];
        double te = tl[i];
        double ti = te - tb;  
        solve (tb,te);
        adapt();
        tb = te;
      }
    }  
  } 
  flush_history();
  *lsb << "done" << endl;
}


void ALGO_WAVE::set_nadapt (sizet nadapt) {
  _nadapt = nadapt;
}

void ALGO_WAVE::set_hist_file_prefix (string hist_file_prefix) {
  if (_hist_file_prefix == "" &&
      hist_file_prefix != "")
    _hist_file_prefix = hist_file_prefix;
}



void  
ALGO_WAVE::output() 
{
  //write history data 
  if ((!_is_back_tracking) || (_force_back_track_output) ) {
    if (_time_between_outputs != 0) {
      if (_t >= _time_between_outputs*_n_output_times ) {
        update_history();
        _n_output_times ++;      
      } 
    } else {
      update_history();  
    }
    // write xhexvue data
    if (_time_between_graphics != 0) {
      if (_t >= _time_between_graphics*_n_graphics_times ) {
        ALGO_XHEXVUE a(this->name(), this->mgr());
        a.make_file_name(_n_graphics_times);
        a.write(_proto_wave,_t);
        _n_graphics_times ++;      
      } 
    } else {
      ALGO_XHEXVUE a(this->name(), this->mgr());
      a.make_file_name(_ncycle);
      a.write(_proto_wave); 
      _ncycle++;
    }
  } 

  if (_algo_motion_pict) {
    update_geometry ();
    _algo_motion_pict->update(_updated_geometry);
  }
}

void  
ALGO_WAVE::increment_t() 
{
  if (_algo_motion_pict) {
    char buf[1024];
    sprintf(buf, "ALGO_WAVE::increment_t(); time=%g", _t);
    _algo_motion_pict->event_loop(false, buf);
  }
  if (_accumulate_error) assemble_error();
  _t += _dt; 
}

void 
ALGO_WAVE::update_history() 
{
  for (list<WATCH_POINT_ECELL*>::iterator i= _watch_point.begin(); i != _watch_point.end(); i++) {
    ECELL *ecell = (*i)->ecell();
    POINT param_loc = (*i)->param_loc();
    ecell->set_watch_point_value ((*i), param_loc);
  }
  double k_energy =  _proto_wave->total_kinetic_energy(); double s_energy = _proto_wave->total_strain_energy();
  FIXED_VECTOR<3> m =  _proto_wave->momentum();
  double mass = _proto_wave->total_mass();
  _global_watch_point->set_time_value("kinetic_energy",k_energy,_t);
  _global_watch_point->set_time_value("strain_energy",s_energy,_t);
   _global_watch_point->set_time_value("total_energy",k_energy + s_energy,_t);
  _global_watch_point->set_time_value("momentum_x",m(0) ,_t);
  _global_watch_point->set_time_value("momentum_y",m(1) ,_t);
  _global_watch_point->set_time_value("momentum_z",m(2) ,_t);
  _global_watch_point->set_time_value("nbfuns",_u->npairs(),_t);
  _global_watch_point->set_time_value("mass",mass,_t);
}

void 
ALGO_WAVE::flush_history() 
{
  for (list<WATCH_POINT_ECELL*>::iterator i= _watch_point.begin(); i != _watch_point.end(); i++) {
    (*i)->flush (); 
  }
}

ALGO_WAVE::~ALGO_WAVE () {
  if (_proto_wave) delete _proto_wave;
  _gsubmeshes.clear();
  if (_algo_refine) delete _algo_refine;
  if (_algo_errest) delete _algo_errest;
  if (_geometry) delete _geometry;
  if (_u) delete _u;
  if (_v) delete _v;
  //if (_ebc) delete _ebc; //RAW: ebc is a smart_handle
  for (list<WATCH_POINT_ECELL*>::iterator i= _watch_point.begin(); i != _watch_point.end(); i++) {
    delete *i;
  }
  _watch_point.clear();
  if (_u0) delete _u0;
  if (_v0) delete _v0;
  if (_updated_geometry) delete _updated_geometry;
  if (_algo_motion_pict) {
    _algo_motion_pict->clear();
    delete _algo_motion_pict;
  }
}

