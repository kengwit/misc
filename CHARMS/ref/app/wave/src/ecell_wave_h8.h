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
#ifndef ECELL_WAVE_H8_H
# define ECELL_WAVE_H8_H

#include "ecell.h"
#include "ecell_wave.h"
#include "field.h"
#include "proto_wave.h"
#include "gaussP.h"
#include "evalpt.h"
#include "func.h"
#include "watchpoint_ecell.h"

class ECELL_WAVE_H8 : public ECELL_WAVE {

 public: // class functions /////////////////////////////////////////

  static bool register_make_func ();

 public: // object functions ////////////////////////////////////////

  ECELL_WAVE_H8 (GCELL *gcell);

  double suggested_time_step(); 

  bool assem_internal_forces ();
  bool assem_body_load ();
  bool assem_lumped_mass();
  double strain_energy (); 
  // void set_watch_point_value (WATCH_POINT *watch_point) ;   
  bool assem_consistent_mass() { return false; }
  void set_watch_point_value (WATCH_POINT *watch_point, POINT &param_loc) ;   
  bool assem_tractions_load () {return true;}
  FIXED_VECTOR<3> calc_loc (POINT &param_loc);
  vector<double> calc_quants (POINT &param_loc,list<string>vars);
  bool calc_stress (POINT &param_loc, double stress[6],double strain[6]);
  void attach (PROTO_WAVE *proto, string ggpath){ ECELL_WAVE::attach (proto,ggpath);}
  #if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION 
  FIXED_VECTOR<3> momentum(){}; 
  double kinetic_energy() {return 0;}
#endif
  #if defined(IMPLICIT_SOLVER) && IMPLICIT_SOLVER
  bool assem_m_bt_k (double m_bt_k[3][3]){return 0;}
  bool assem_bt_k(double K[3][3] ){return 0;}
  bool assem_k_d_bar (){return 0;}
#endif
 private: // class data /////////////////////////////////////////////

  static const char _type_name[];
 
};

#endif
