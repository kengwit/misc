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
#ifndef ECELL_ERREST_WAVE_H8_H
# define ECELL_ERREST_WAVE_H8_H

#include "ecell.h"
#include "ecell_errest.h"
#include "field.h"
extern "C" {
#include "gaussP.h"
}
#include "evalpt.h"
#include "mat_viscel.h"
#include "mat_viscel_iso.h"
#include "watchpoint.h"

class ECELL_ERREST_WAVE_H8 : public ECELL_ERREST<3> {

 public: // class functions /////////////////////////////////////////

  static bool register_make_func ();

 public: // object functions ////////////////////////////////////////

  /**
   */
  ECELL_ERREST_WAVE_H8 (GCELL *gcell);
  
  /**
   */
  void assemble_error (FIELD_VECTOR *u, FIELD_VECTOR *geometry, class PROTO_ERREST<3> *proto);
  void attach (class PROTO_ERREST<3> *proto, string ggpath);
  /**
     p = Int(Rho*V*d_omega)
       = Int(Rho*Ni*Vi*d_omega)
       = Vi*Int(Rho*Ni*d_omega)
       = Vi*Mi
  */

  double strain_energy (){ return 0; }; 


 private: // class data /////////////////////////////////////////////
  
  static const char TYPE_NAME[];
  MAT_VISCEL *_mat;

};

#endif
