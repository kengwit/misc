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
#ifndef ECELL_WAVE_TRACTION_L2_H
# define ECELL_WAVE_TRACTION_L2_H

#include "ecell.h"
#include "ecell_wave.h"
#include "field.h"
#include "gaussP.h"
#include "evalpt.h"
#include "func.h"
#include "watchpoint_ecell.h"
#include "proto_wave.h"

class ECELL_WAVE_TRACTION_L2 : public ECELL_WAVE {

 public: // class functions /////////////////////////////////////////

  static bool register_make_func ();

 public: // object functions ////////////////////////////////////////

  ECELL_WAVE_TRACTION_L2 (GCELL *gcell);
  
  bool assem_internal_forces () {return true;}
  bool assem_body_load () {return true;}
  bool assem_lumped_mass() {return true;}
  bool assem_consistent_mass() { return true; }
  bool assem_tractions_load ();
  void attach (PROTO_WAVE *proto, string ggpath){ ECELL_WAVE::attach (proto,ggpath);}
  double strain_energy (){return 0;} 
  #if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION 
  FIXED_VECTOR<3> momentum(){return 0;}
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
