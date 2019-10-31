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
#include "proto_wave.h"
#include "mat_viscel.h"
#include "algo_wave.h"
#include "ecell_wave_h8.h"
#include "select_mass.h"

using namespace std;
// How to ensure that an algorithm uses the protocol for
// a field defined over the proper set of submeshes, and of the
// proper meaning (scalar field)?
PROTO_WAVE::PROTO_WAVE (ALGO_WAVE *algo_wave,// algorithm using this protocol
                        DB *db, // database
                        FIELD_VECTOR *geometry,
                        FIELD_VECTOR *u,
                        FIELD_VECTOR *v
                        ) 
{
  // instantiated for (with)
  _db = db; CHECK (_db != NULL, EXCEPTION_NULL_PTR,;);
  _mat_map = new MAT_MAP <MAT_VISCEL> (algo_wave->mgr()->mat_mgr());
  _algo_wave = algo_wave;
  CHECK (_algo_wave != NULL, EXCEPTION_NULL_PTR,;);
  _u = u; _v = v;
  CHECK (_u != NULL, EXCEPTION_NULL_PTR,;);
  CHECK (_v != NULL, EXCEPTION_NULL_PTR,;);
  CHECK(_u->bfun_set()==_v->bfun_set(), EXCEPTION_ILLEGAL_USE,;);
  _geometry = geometry;
  CHECK (_geometry != NULL, EXCEPTION_NULL_PTR,;);
  // now create the ecells
  _ecells.clear();
  BFUN_SET::gsubmesh_enumerator_t sme = _u->bfun_set()->gsubmesh_enumerator ();
  GSUBMESH *gsubmesh;
  sme.reset ();
  int temp_counter = 0;
  while ((gsubmesh = sme.next ())) {
    // gsubmesh->debug_display ();
    GSUBMESH::gcell_group_enumerator_t gge = gsubmesh->gcell_group_enumerator ();
    GCELL_GROUP *gcell_group;
    gge.reset ();
    while ((gcell_group = gge.next ())) {
      string ggpath = "algorithms/wave/"
        + algo_wave->name () + "/gcell_groups/" + gcell_group->name ();
      string implementation = db->DB_GET_STRING (ggpath + "/implementation");
      GCELL_GROUP::gcell_enumerator_t ge = gcell_group->gcell_enumerator ();
      GCELL *gcell;
      ge.reset ();  
      while ((gcell = ge.next ())) {
        if (gcell->is_topmost_gcell_in_field (_u)) {
          ECELL_WAVE *ecell = ECELL_WAVE::make_ecell (gcell, implementation); 
          if (!ecell) {
            cerr<< "failed to create ecell\t"<< gcell->type_name() << "\t"  <<implementation<<"\n";
            CHECK (ecell, EXCEPTION_NULL_PTR,;);
          }
          ecell->attach (this, ggpath);
          _ecells.push_back (ecell);  
          // cerr << "Made ecell for "; gcell->debug_display (); cerr << endl;
        }
      }
    }
  }
  _fm = 0;
}



bool PROTO_WAVE::assemble_eff_loads ()
{
  for (sizet j=0; j<_u->npairs(); j++) {
    _fm[j].force_t = _fm[j].force_t_dt;
    _fm[j].force_t_dt = 0;
  }
  for (list<ECELL_WAVE *>::const_iterator ecell = _ecells.begin ();
       ecell != _ecells.end (); ecell++) {
    if (!(*ecell)->assem_body_load ()) return false;
    if (!(*ecell)->assem_tractions_load ()) return false;
    if (!(*ecell)->assem_internal_forces ()) return false;
   }


#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION
  _slesctx->zero_rhs (_slesctx);
  for (sizet j=0; j<_u->npairs(); j++) {
    double v[3] = { _fm[j].force_t_dt(0), _fm[j].force_t_dt(1), _fm[j].force_t_dt(2) };
    int rows[3] = { j*3+0, j*3+1, j*3+2 };
    _slesctx->assemble_vect (_slesctx, FLES_ADD, 3, rows, v);
  }
  _slesctx->finish_rhs (_slesctx);
  //  _slesctx->dump_lhs (_slesctx, "lhsdump_eff"); //RAW
  if (!_slesctx->solve (_slesctx)) {
    // RAW: recovery?
    CHECK(0, EXCEPTION_NULL_PTR,;);
  }
  for (sizet j=0; j<_u->npairs(); j++) {
    for (sizet c=0; c < 3; c++) {
      _fm[j].force_t_dt(c) = _slesctx->get_from_sol (_slesctx, j*3+c);
    }
  }
#endif
  
  return true;
}

#if defined(IMPLICIT_SOLVER) && IMPLICIT_SOLVER
bool PROTO_WAVE::calc_load()
{
  for (list<ECELL_WAVE *>::const_iterator ecell = _ecells.begin ();
       ecell != _ecells.end (); ecell++) {
    if (!(*ecell)->assem_body_load ()) return false;
    if (!(*ecell)->assem_tractions_load ()) return false;
    if (!(*ecell)->assem_internal_forces ()) return false;
   }
}


bool PROTO_WAVE::implicite_begin_load_calc ()
{
  calc_load();
#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION
  _slesctx->zero_rhs (_slesctx);
  for (sizet j=0; j<_u->npairs(); j++) {
    double v[3] = { _fm[j].force_t_dt(0), _fm[j].force_t_dt(1), _fm[j].force_t_dt(2) }; //force_t_dt stores A0
    int rows[3] = { j*3+0, j*3+1, j*3+2 };
    _slesctx->assemble_vect (_slesctx, FLES_ADD, 3, rows, v);
  }
  _slesctx->finish_rhs (_slesctx);
  //  _slesctx->dump_lhs (_slesctx, "lhsdump_eff"); //RAW
  if (!_slesctx->solve (_slesctx)) {
    // RAW: recovery?
    CHECK(0, EXCEPTION_NULL_PTR,;);
  }
  for (sizet j=0; j<_u->npairs(); j++) {
    for (sizet c=0; c < 3; c++) {
      _fm[j].force_t_dt(c) = _slesctx->get_from_sol (_slesctx, j*3+c);
    }
  }
#endif
  return true;
}
#endif

#if defined(IMPLICIT_SOLVER) && IMPLICIT_SOLVER
bool PROTO_WAVE::implicite_assemble_eff_loads ()
{
  for (sizet j=0; j<_u->npairs(); j++) {
    _fm[j].force_t = _fm[j].force_t_dt; //transfer A0 to force_t
    _fm[j].force_t_dt = 0;
  }
  // equivalent to F-K*D_tilda; F than d_bar and v_var and Kd_bar and An+1

 for (sizet j=0; j<_u->npairs(); j++) { //RAW this is not efficient implementation. here accelarations can be directly stored
    _force[j].force_t_dt = _fm[j].force_t_dt;
  }

  //think about when it comes here in next round
  for (list<ECELL_WAVE *>::const_iterator ecell = _ecells.begin ();
       ecell != _ecells.end (); ecell++) {  //push force in f_t_dt
    if (!(*ecell)->assem_k_d_bar ()) return false;
   }
  _slesctx1->zero_rhs (_slesctx);
  for (sizet j=0; j<_u->npairs(); j++) {
    double v[3] = { _fm[j].force_t_dt(0), _fm[j].force_t_dt(1), _fm[j].force_t_dt(2) };
    int rows[3] = { j*3+0, j*3+1, j*3+2 };
    _slesctx1->assemble_vect (_slesctx1, FLES_ADD, 3, rows, v);
  }
  _slesctx1->finish_rhs (_slesctx1);
  //  _slesctx->dump_lhs (_slesctx, "lhsdump_eff"); //RAW
  if (!_slesctx1->solve (_slesctx1)) {
    // RAW: recovery?
    CHECK(0, EXCEPTION_NULL_PTR,;);
  }
  for (sizet j=0; j<_u->npairs(); j++) {
    for (sizet c=0; c < 3; c++) {
      _fm[j].force_t_dt(c) = _slesctx1->get_from_sol (_slesctx1, j*3+c);
    }
  }
  // An+1 is in force_t_dt
  return true;
}
#endif

#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION
static int _tot_of_equations (void *client)
{
  PROTO_WAVE *p = (PROTO_WAVE *)client;
  int n = p->u()->npairs()*3; 
  return n;
}

static int _error_notify (fles_error_enum error_num,
                          int line,
                          char *file,
                          char *detail,
                          void *client_data)
{
  char buf[1024];
  sprintf (buf, "%s (l.%d), %s", file, line, detail);
  CHECK (0, EXCEPTION_IN_LOCAL_LIBS, (buf));
  return 0;
}
#endif


bool PROTO_WAVE::assemble_mass ()
{
#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION
  SMART_HANDLE<LOGGER_STREAM> lsb = this->_algo_wave->mgr()->logger("assemble consistent mass", true);
  {
    double nz_per_row = 0.01;
    _slesctx = fles_new_spetsc_solver (this, _tot_of_equations, _error_notify);
  }

  _slesctx->zero_lhs (_slesctx);
  for (list<ECELL_WAVE *>::const_iterator ecell = _ecells.begin ();
       ecell != _ecells.end (); ecell++) {
    if (!(*ecell)->assem_consistent_mass()) return false;
  }
  _slesctx->finish_lhs (_slesctx, FLES_FINAL_FINISH);
  //  _slesctx->dump_lhs (_slesctx, "lhsdump_ass"); //RAW 
#else
  SMART_HANDLE<LOGGER_STREAM> lsb = this->_algo_wave->mgr()->logger("assemble lumped mass", true);
  for (sizet j=0; j<_u->npairs(); j++) {
    _fm[j].mass = 0;
  }
  for (list<ECELL_WAVE *>::const_iterator ecell = _ecells.begin ();
       ecell != _ecells.end (); ecell++) {
    if (!(*ecell)->assem_lumped_mass()) return false;
  }
#endif
  *lsb << "done"  << endl;
  return true;
}

#if defined(IMPLICIT_SOLVER) && IMPLICIT_SOLVER
bool PROTO_WAVE::assem_m_bt_k() {
  _slesctx1 = fles_new_spetsc_solver (this, _tot_of_equations, _error_notify);
  _slesctx1->zero_lhs (_slesctx1);
 
  for (list<ECELL_WAVE *>::const_iterator ecell = _ecells.begin ();
       ecell != _ecells.end (); ecell++) {
    double m_bt_k[3][3]   = {{0,0,0},{0,0,0},{0,0,0}}; //double 3x3
    if (!(*ecell)->assem_bt_k(m_bt_k)) return false;
    if (!(*ecell)->assem_m_bt_k (m_bt_k)) return false;
  }
  _slesctx1->finish_lhs (_slesctx1, FLES_FINAL_FINISH);
}
#endif

void PROTO_WAVE::init_fm () {
  _fm = new fm[_u->npairs()];
  for (sizet j=0; j<_u->npairs();j++) {
   _fm[j].force_t = 0;
   _fm[j].force_t_dt = 0;
#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION
#else
   _fm[j].mass    = 0;
#endif
  }
}

#if defined(IMPLICIT_SOLVER) && IMPLICIT_SOLVER
void PROTO_WAVE::implicite_init_fm () {
  _fm = new fm[_u->npairs()];
  for (sizet j=0; j<_u->npairs();j++) {
   _fm[j].force_t = 0;
   _fm[j].force_t_dt = 0;
#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION
    _force[j].force_t_dt = 0;
#else
   _fm[j].mass    = 0;
#endif
  }
}  
#endif

bool PROTO_WAVE::advance (sizet nsteps) {
  init_fm();
  if(!assemble_mass()) return false;
  SMART_HANDLE<LOGGER_STREAM> lsb = this->_algo_wave->mgr()->logger("advance", true);
  if(!assemble_eff_loads()) return false;
  for (sizet i = 0; i< nsteps; i++){
    _algo_wave->output();
    update_u ();
    assemble_eff_loads();
    update_v();
    _algo_wave->increment_t();
  }  
  _algo_wave->output();
  *lsb << "done"  << endl;
  return true;
}


#if defined(IMPLICIT_SOLVER) && IMPLICIT_SOLVER
bool PROTO_WAVE::implicite_advance (sizet nsteps) {
  implicite_init_fm();
  if(!assemble_mass()) return false;
  SMART_HANDLE<LOGGER_STREAM> lsb = this->_algo_wave->mgr()->logger("advance", true);
  
  if(!implicite_begin_load_calc ()) return false; //fm.force_t_dt  has A0 now 
  for (sizet i = 0; i< nsteps; i++){
    _algo_wave->output();
    assem_m_bt_k(); // equivalent to M+beta*dt**2*k
    implicite_assemble_eff_loads(); // equivalent to F-K*D_tilda; F than d_bar and v_var and Kd_bar and An+1 
    implicite_update_u ();
    update_v();
    _algo_wave->increment_t();
    calc_load();
  }  
  _algo_wave->output();
  *lsb << "done"  << endl;
  return true;
}
#endif

void PROTO_WAVE::print_energy (){
  double str_energy =  total_strain_energy(); double k_energy = total_kinetic_energy();
  SMART_HANDLE<LOGGER_STREAM> lsb = this->_algo_wave->mgr()->logger("print_energy", true);
  *lsb <<"************************************************"<<endl_notime;
  *lsb << "total mass = " << total_mass() << endl_notime;
  *lsb << "total kinetic energy =\t " << k_energy << endl_notime;
  *lsb << "total_momentum \t" << moment() << endl_notime;
  *lsb << "total_strain_energy \t" << str_energy << endl_notime;
  *lsb << "total_energy \t" << (str_energy + k_energy) << endl_notime;
  *lsb <<"************************************************"<<endl_notime;
}


bool PROTO_WAVE::update_u(){
  double dt = _algo_wave->dt();
  double dt2 = dt*dt/2.0;
  // debug_display();
  for (sizet j=0; j<_u->npairs();j++) {
    FIELD_PAIR<3> *fp  = _u->field_pair(j);
    FIXED_VECTOR<3> ut = fp->dofparam();
    FIXED_VECTOR<3> vt = _v->field_pair(j)->dofparam();
#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION
    ut.add (dt, vt, dt2, _fm[j].force_t_dt);
#else
    ut.add (dt, vt, dt2/_fm[j].mass, _fm[j].force_t_dt);
#endif
    fp->set_dofparam(ut);
    //  cout<<"ut=\t"<<ut<<"\n";
  }
  _algo_wave->ebc()->apply(_u, _algo_wave->t());
  return true;
} 

#if defined(IMPLICIT_SOLVER) && IMPLICIT_SOLVER
bool PROTO_WAVE::implicite_update_u(){
  double dt = _algo_wave->dt();
  double dt4 = dt*dt/4.0;
  for (sizet j=0; j<_u->npairs();j++) {
    FIELD_PAIR<3> *fp  = _u->field_pair(j);
    FIXED_VECTOR<3> ut = fp->dofparam();
    FIXED_VECTOR<3> vt = _v->field_pair(j)->dofparam();
    ut.add (dt, vt, dt4, (_fm[j].force_t_dt +_fm[j].force_t) );
    fp->set_dofparam(ut);
  }
  return true;
}
#endif

bool PROTO_WAVE::update_v (){
  double dt = _algo_wave->dt();
  //_v->enforce_ebc();
  for (sizet j=0; j<_v->npairs();j++) {
    FIELD_PAIR<3> *fp  = _v->field_pair(j);
    FIXED_VECTOR<3> vt = fp->dofparam();
#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION
    double dm = dt/(2);
    vt.add (dm, _fm[j].force_t_dt, dm,_fm[j].force_t); 
#else
    double dm = dt/(2*_fm[j].mass);
    vt.add (dm, _fm[j].force_t_dt, dm,_fm[j].force_t); 
#endif
    fp->set_dofparam(vt);
  }
  return true;
} 


double PROTO_WAVE::stable_time_step()
{
  double ts,ts_temp ; ts = 0; ts_temp = FLT_MAX;
  for (list<ECELL_WAVE *>::iterator i = _ecells.begin (); i != _ecells.end (); i++) {
    if ((*i)) ts = (*i)->suggested_time_step();
    if (ts < ts_temp) ts_temp = ts;   
  }
  return ts_temp;
}



MAT_VISCEL *PROTO_WAVE::mat (string mattype, string matname)
{
  return _mat_map->mat (mattype, matname);
}

LOAD_ON_GCELL *PROTO_WAVE::load (string loadtype, string loadname)
{
  return dynamic_cast <LOAD_ON_GCELL *>(this->_algo_wave->mgr()->load_mgr()->load (loadtype, loadname));
}

PROTO_WAVE::~PROTO_WAVE () {
  for (list<ECELL_WAVE  *>::const_iterator i = _ecells.begin ();
       i != _ecells.end (); i++) {
    ECELL_WAVE *ecell = *i;
    delete ecell;
  }
  _ecells.clear();
  if (_mat_map) delete _mat_map;
  if (_fm) delete [] _fm;
#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION
  if (_slesctx) _slesctx->destroy_ctx (_slesctx);
  _slesctx = 0;
#endif
}  


void PROTO_WAVE::assemble_m(BFUN_DOFPARAM_PAIR_ID dpid,double value )
{
#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION
#else
  CHECK(!isnan(value), EXCEPTION_BAD_VALUE,;);
  _fm[dpid].mass += value; 
#endif
}


void PROTO_WAVE::assemble_f(BFUN_DOFPARAM_PAIR_ID dpid,double value[3] )
{
 for (sizet i=0; i<3; i++)_fm[dpid].force_t_dt(i) += value[i];
}

#if defined(IMPLICIT_SOLVER) && IMPLICIT_SOLVER
void PROTO_WAVE::assemble_f_kd(BFUN_DOFPARAM_PAIR_ID dpid,double value[3] )
{
 for (sizet i=0; i<3; i++)_fm[dpid].force_t_dt(i) += value[i];
}
#endif

#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION
void PROTO_WAVE::assemble_m (BFUN_DOFPARAM_PAIR_ID dpid1, BFUN_DOFPARAM_PAIR_ID dpid2, double value)
{
  double m[3][3]={{0,0,0},{0,0,0},{0,0,0}};
  m[0][0] = m[1][1] = m[2][2] = value;
  int rows[3] = { dpid1*3+0, dpid1*3+1, dpid1*3+2 }; 
  int cols[3] = { dpid2*3+0, dpid2*3+1, dpid2*3+2 }; 
  _slesctx->assemble_mat (_slesctx, FLES_FALSE, FLES_ADD, 3, rows, 3, cols, (double *)m);
}
#endif

#if defined(IMPLICIT_SOLVER) && IMPLICIT_SOLVER
void PROTO_WAVE::assemble_m_bt_k (BFUN_DOFPARAM_PAIR_ID dpid1, BFUN_DOFPARAM_PAIR_ID dpid2, double value, double m[3][3])
{
  m[0][0]  += value; 
  m[1][1]  += value;
  m[2][2]  += value;
  int rows[3] = { dpid1*3+0, dpid1*3+1, dpid1*3+2 }; 
  int cols[3] = { dpid2*3+0, dpid2*3+1, dpid2*3+2 }; 
  _slesctx1->assemble_mat (_slesctx, FLES_FALSE, FLES_ADD, 3, rows, 3, cols, (double *)m);
}
#endif


double  PROTO_WAVE::time() {return (_algo_wave->t());} 

void PROTO_WAVE::debug_display() {
#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION
#else
  cout << "force1"<<"\t"<< "force2"<<"\t"<< "force3"<<"\t"<< "mass"<<"\n";
  for (sizet i=0; i < _u->npairs(); i++) {
    cout <<_fm[i].force_t_dt(0)<<"\t"<<_fm[i].force_t_dt(1)<<"\t"<<_fm[i].force_t_dt(2)<<"\t"<<_fm[i].mass<<"\n";
  }
#endif
}

void 
PROTO_WAVE::map_watchpoints (list<WATCH_POINT_ECELL*> watch_point) 
{
  for (list<WATCH_POINT_ECELL*>::iterator i= watch_point.begin(); i != watch_point.end(); i++) {
     for (list<ECELL_WAVE *>::const_iterator ecell = _ecells.begin ();
          ecell != _ecells.end () ; ecell++) {
       if ( (*ecell)->gcell()->root()->id() == (*i)->gcell_id()) {
         if ((*i)->set_ecell(*ecell)) goto mapped; 
       } 
     }
            CHECK(0, EXCEPTION_NULL_PTR,;); // watchpoint could not be mapped!? 
  mapped:;
  }
}

#if defined(IMPLICIT_SOLVER) && IMPLICIT_SOLVER
double PROTO_WAVE::dt() {return _algo_wave->dt();}
#endif

void PROTO_WAVE::add_graphics (ALGO_MOTION_PICT *algo_motion_pict, FIELD_VECTOR *geometry) {
  for (list<ECELL_WAVE  *>::const_iterator i = _ecells.begin ();
       i != _ecells.end (); i++) {
    ECELL_WAVE *ecell = *i;
    algo_motion_pict->add_graphics(ecell->gcell(), geometry);
  }
  for (sizet j=0; j<_u->npairs();j++) {
    algo_motion_pict->add_graphics(geometry->field_pair(j)->bfun()->fen(), geometry);
  }  
}
