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
#ifndef PROTO_WAVE_H
# define PROTO_WAVE_H

#include <list>
#include "algo.h"
#include "gmesh.h"
#include "ecell_wave.h"
#include "mat_mgr.h"
#include "mat_map.h"
#include "mat_viscel.h"
#include "load_on_gcell.h"
#include "watchpoint_ecell.h"
#include "quant_evaluator.h"
#include "proto_base.h"
#include "select_mass.h"
#include "algo_motion_pict.h"
#include "select_solver.h"

/**
  Prototypical wave equation: M u" + K u = F

  K=stiffness matrix
  u=displacement vector
  F=load vector

 */
class PROTO_WAVE : public PROTO_BASE {

 public: // object functions ////////////////////////////////////////

   PROTO_WAVE (class ALGO_WAVE *algo_wave,
               DB *db, FIELD_VECTOR *geometry,
               FIELD_VECTOR *u,FIELD_VECTOR *v);
  ~PROTO_WAVE () ;
  /**
     Assemble the stiffness matrix K (see the class comment).
  */
  bool assemble_restoring_force ();
  /**
     Assemble the load (rhs) terms corresponding to non-zero prescribed
     boundary displacements, applied body loads, and applied tractions.
  */
  bool assemble_load_terms ();
  /**
    
  */
  bool advance(sizet  nsteps); 
                
  void add_graphics (ALGO_MOTION_PICT *algo_motion_pict, FIELD_VECTOR *geometry);

  double stable_time_step();
  /**
     Get the primary field that this protocol is used
     to solve for.This will be the latest displacement.
  */
  //  FIELD_VECTOR *u () { return _u; }//RAW
  /**
     Get the geometry field that this protocol is supposed
     to use in the solution.
  */
  // FIELD_VECTOR *geometry () { return _geometry; } //RAW
  /**
     Get the associated database.
  */
  DB *db () { return _db; }
  /**
     Material.
  */
  MAT_VISCEL *mat (string mattype, string matname) ;
  /**
     Load.
  */
  LOAD_ON_GCELL *load (string loadtype, string loadname) ;
  
  
  void assemble_f (BFUN_DOFPARAM_PAIR_ID dpid, double value[3]);
  
  void assemble_m (BFUN_DOFPARAM_PAIR_ID dpid, double value);

  void assemble_m (BFUN_DOFPARAM_PAIR_ID dpid1, BFUN_DOFPARAM_PAIR_ID dpid2, double value);

  bool assemble_eff_loads ();

  bool assemble_mass ();
  void print_energy();
  void init_fm ();
  bool update_u();
  bool update_v();
  double time();
  void   debug_display();
  void   map_watchpoints(list<WATCH_POINT_ECELL*> watch_point);
  FIELD_VECTOR *u() { return _u;}
  FIELD_VECTOR *v() { return _v;}
  FIELD_VECTOR *geometry() {return _geometry;}
 public:
#if defined(IMPLICIT_SOLVER) && IMPLICIT_SOLVER
  FIXED_VECTOR<3> get_acc_t    (sizet j) { return _fm[j].force_t;}
  FIXED_VECTOR<3> get_acc_t_dt (sizet j) { return _fm[j].force_t_dt;}
  void assemble_m_bt_k (BFUN_DOFPARAM_PAIR_ID dpid1, BFUN_DOFPARAM_PAIR_ID dpid2, double value, double m[3][3]);
  void assemble_f_kd(BFUN_DOFPARAM_PAIR_ID dpid,double value[3] );
  bool implicite_update_u();
  bool implicite_advance (sizet nsteps);
  bool assem_m_bt_k();
  bool implicite_assemble_eff_loads ();
  bool implicite_begin_load_calc ();
  bool calc_load();
  void implicite_init_fm () ;
  double dt();
#endif 

 public:

  class QUANT_EVAL_WAVE : public QUANT_EVALUATOR {
    
  public:
    
    QUANT_EVAL_WAVE (ECELL  *e, list<string> vars, PROTO_WAVE *p): QUANT_EVALUATOR (e, vars) {
      _e = dynamic_cast <ECELL_WAVE *> (e);
      _p = p;
      for (list<string>::iterator i = vars.begin(); i != vars.end(); i++) _vars.push_back(*i); 
    }
    
    FIXED_VECTOR<3> loc (POINT &param_loc) {
      return  _e->calc_loc (param_loc);
      
    }
    vector <double> quants (POINT &param_loc) {
      return _e->calc_quants(param_loc, _vars); 
    }
    
    private :

      ECELL_WAVE     *_e;
      list<string>    _vars;
      PROTO_WAVE     *_p; 

  };
 
  SMART_HANDLE<QUANT_EVALUATOR> quant_evaluator (ECELL *e,list<string>  quants) {
    SMART_HANDLE<QUANT_EVAL_WAVE> q = new QUANT_EVAL_WAVE(e, quants, this);
    return q;    
  }

  
  class ECELL_ENUMERATOR_WAVE : public ECELL_ENUMERATOR  {
  public:
    ECELL_ENUMERATOR_WAVE (PROTO_WAVE *p): ECELL_ENUMERATOR(p) { 
      _e = new ENUMERATOR <list <ECELL_WAVE *>, ECELL_WAVE>(p->ecells()) ;
    } 
    virtual ~ECELL_ENUMERATOR_WAVE() {
      delete (_e);
     };
    ECELL *next () {
      return (_e->next());     
      }
    void reset() {
       return (_e->reset());
    } 
  private:
    ENUMERATOR <list <ECELL_WAVE *>, ECELL_WAVE> *_e;
  };
 

   ECELL_ENUMERATOR *ecell_enumerator () {
     ECELL_ENUMERATOR_WAVE *e =  new ECELL_ENUMERATOR_WAVE(this);
     return e;
   };

   list <ECELL_WAVE *> ecells() {return _ecells;}

 private: // force buffer type
  
  typedef struct {
    FIXED_VECTOR<3>               force_t;
    FIXED_VECTOR<3>               force_t_dt;
#if defined(CONSISTENT_MASS_ADAPTION) && (!CONSISTENT_MASS_ADAPTION)
    double                        mass;
#endif
  } fm;

#if defined(IMPLICIT_SOLVER) && IMPLICIT_SOLVER
 typedef struct {
    FIXED_VECTOR<3>               force_t_dt;
 } force;
#endif
  
 private: // object data ///////////////////////////////////////////// 

  class ALGO_WAVE                *_algo_wave;
  DB                             *_db;
  list < ECELL_WAVE *>            _ecells; 
  MAT_MAP <MAT_VISCEL >          *_mat_map;
  fm                             *_fm;
#if defined(IMPLICIT_SOLVER) && IMPLICIT_SOLVER
  force                          *_force;
#endif
  FIELD_VECTOR                   *_u;
  FIELD_VECTOR                   *_v;
  FIELD_VECTOR                   *_geometry;
#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION
  fles_solver_ctx_t              _slesctx;
  fles_solver_ctx_t              _slesctx1;
#endif

 public:

  double total_mass () { // RAW this is not correct: total mass = sum_I rho_I integral(N_I)
    double tot_mass = 0;
#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION
#else
    for (sizet j=0; j<_u->npairs(); j++) {
      tot_mass += _fm[j].mass;
    }
#endif
    return tot_mass;
  }
  

double total_kinetic_energy () {
    double K = 0;
#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION
   for (list<ECELL_WAVE*>::iterator i = _ecells.begin(); i != _ecells.end(); i++) {
    K += (*i)->kinetic_energy();
 }  
#else
    for (sizet j=0; j<_u->npairs(); j++) {
      FIXED_VECTOR<3> vt = _v->jth_field_pair(j)->dofparam();
      K += _fm[j].mass * vt.l2_norm_squared();
    }
    K = K/2;
#endif
    return K;
  }

FIXED_VECTOR<3> momentum () {
 FIXED_VECTOR<3> p(0);

#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION
 for (list<ECELL_WAVE*>::iterator i = _ecells.begin(); i != _ecells.end(); i++) {
   p += (*i)->momentum();
 }  
#else
     for (sizet j=0; j<_u->npairs(); j++) {
       FIXED_VECTOR<3> vt = _v->jth_field_pair(j)->dofparam();
       // p.sadd( _fm[j].mass, vt);    //RAW 
       p.add( _fm[j].mass, vt); 
     }
#endif
     return p;
  }


double moment () {
#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION
#else
   double Pz = 0;
    for (sizet j=0; j<_u->npairs(); j++) {
      FIXED_VECTOR<3> vt = _v->jth_field_pair(j)->dofparam();
      Pz +=  _fm[j].mass * (vt(0)+vt(1)+vt(2));
    }
    // FIXED_VECTOR <3> Mz  = momentum();
    // if (Pz != (Mz(0)+ Mz(1) + Mz(2))) cout<<"momentums are not equal. Mz = "<<  (Mz(0)+ Mz(1) + Mz(2)) << "\n"; 
   return Pz;
#endif
    return 0;
 }

  double total_strain_energy() {
    double strain_energy = 0;
    for (list <ECELL_WAVE *>:: iterator i = _ecells.begin(); i != _ecells.end(); i++) {
      strain_energy += (*i)->strain_energy();
    }  
    return strain_energy;
  }  
};

#endif


/*
 double momentum_x () {
#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION
#else
   double Px = 0;
    for (sizet j=0; j<_v->npairs(); j++) {
      FIXED_VECTOR<3> vt = _v->jth_field_pair(j)->dofparam();
      Px+ =  _fm[j].mass * vt(0);
    }
#endif
    return Px;
 } 

 double momentum_y () {
#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION
#else
   double Py = 0;
    for (sizet j=0; j<_u->npairs(); j++) {
      FIXED_VECTOR<3> vt = _v->jth_field_pair(j)->dofparam();
      Py+ =  _fm[j].mass * vt(1);
    }
#endif
    return Py;
 } 

 double momentum_z () {
#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION
#else
   double Pz = 0;
    for (sizet j=0; j<_u->npairs(); j++) {
      FIXED_VECTOR<3> vt = _v->jth_field_pair(j)->dofparam();
      Pz+ =  _fm[j].mass * vt(z);
    }
#endif
    return Pz;
 }
*/
