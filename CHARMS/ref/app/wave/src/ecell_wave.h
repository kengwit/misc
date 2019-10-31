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
#ifndef ECELL_WAVE_H
# define ECELL_WAVE_H

#include "ecell.h"
#include "field.h"
#include "ofactory.h"
#include "func.h"
#include "mat_viscel.h"
#include "load_on_gcell_body.h"
#include "load_on_gcell_traction.h"
#include "watchpoint_ecell.h"
#include "quant_evaluator.h"
#include "select_mass.h"
#include "select_solver.h"

/**
   The ecells of this class are used by the PROTO_wave
   protocol.
*/
class ECELL_WAVE : public ECELL {

 public: // declarations ////////////////////////////////////////////

  /**
     This is the type of the function that specializations of ECELL_WAVE
     need to register with this class so that they can be manufactured
     by the factory of this class.
  */
  typedef ECELL_WAVE *(*make_ecell_func) (GCELL *gcell);

 public: // class functions  ////////////////////////////////////////

  /**
     Clients of this class may call this function to produce instances
     of the derived classes.  They need to pass an instance of a gcell, and
     a string indicating which implementation of the ecell should be created.
     In other words, given (a) the type of the gcell, and (b) a descriptive name
     (string) that distinguishes among different implementations of the same behaviour,
     the factory of this class will produce an instance of this class (actually,
     a descendant of this class that overrides the required methods).
   */
  static ECELL_WAVE *make_ecell (GCELL *gcell, string implementation);

 protected:
  
  /**
     Call this function to register a make function for the ecell
     based on the gcell type specified, and of given implementation.
  */
  static bool register_make_func (make_ecell_func make, string gcell_type, string implementation);

 public: // object functions ////////////////////////////////////////

  /**
   */
  ECELL_WAVE (GCELL *gcell) : ECELL(gcell) {
    _mat = 0;
    _body_load = 0;
    _traction_load = 0;
  }
  
  virtual bool assem_internal_forces () = 0;
  virtual bool assem_body_load () = 0;
  virtual bool assem_tractions_load () = 0;
  virtual bool assem_lumped_mass ()= 0;
#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION
  virtual bool assem_consistent_mass ()= 0;
#endif
  virtual double suggested_time_step() { return FLT_MAX; }
  virtual FIXED_VECTOR<3> calc_loc (POINT &param_loc) { return FIXED_VECTOR<3> (0); }
  virtual double strain_energy()  { return  0; }
  virtual vector<double> calc_quants (POINT &param_loc, list<string>vars) {
    vector<double> v;
    v.assign(vars.size(), QUANT_EVALUATOR::INVALID_DATA);
    return v;
  }
  virtual bool calc_stress (POINT &param_loc, double stress[6], double strain[6]) { return false; }
  
  /**
     Attach the evaluation cell to the protocol.
     This includes:
     (a) get all needed parameters from the database;
     (b) establish a binding of the degrees of freedom at the nodes
         interacting at the ecell to the LES.
     Arguments: database, les (linear equation solver), field name (what name you wish
     the solver to know the field under), field.
   */
  virtual void attach (class PROTO_WAVE *proto, string ggpath);
  /**
   */
  virtual ~ECELL_WAVE () {};

  /**
     Set the name of material this cell should use.
  */
  void set_mat (MAT_VISCEL *mat) { _mat = mat; }
  /**
     Get the material the ecell should use.
     */
  MAT_VISCEL *mat () const { return _mat; }
  /**
     Get the body load (could be null, in which case no body load is defined).
  */
  LOAD_ON_GCELL_BODY *body_load () const { return _body_load; }
  /**
     Get the traction load (could be null, in which case no traction load is defined).
  */
  LOAD_ON_GCELL_TRACTION *traction_load () const { return _traction_load; }

#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION 
  virtual FIXED_VECTOR<3> momentum() = 0; 
  virtual double kinetic_energy() = 0;
#endif

#if defined(IMPLICIT_SOLVER) && IMPLICIT_SOLVER
  virtual bool assem_m_bt_k (double m_bt_k[3][3]) = 0;
  virtual bool assem_bt_k(double K[3][3] ) = 0;
  virtual bool assem_k_d_bar () = 0;
#endif  
 private: // class objects  /////////////////////////////////////////

  static OFACTORY<ECELL_WAVE, make_ecell_func, GCELL *> *_ecell_factory;

 protected:
  
  class PROTO_WAVE              *_proto;

 private: // object data ////////////////////////////////////////////

  MAT_VISCEL                    *_mat;
  LOAD_ON_GCELL_BODY            *_body_load;
  LOAD_ON_GCELL_TRACTION        *_traction_load;

};

#endif
