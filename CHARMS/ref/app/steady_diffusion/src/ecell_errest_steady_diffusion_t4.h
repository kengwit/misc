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
#ifndef ECELL_ERREST_STEADY_DIFFUSION_T4_H
# define ECELL_ERREST_STEADY_DIFFUSION_T4_H

#include "ecell.h"
#include "ecell_errest.h"
#include "field.h"
extern "C" {
#include "gaussP.h"
}
#include "evalpt.h"

class ECELL_ERREST_STEADY_DIFFUSION_T4 : public ECELL_ERREST<1> {

 public: // class functions /////////////////////////////////////////

  static bool register_make_func ();

 public: // object functions ////////////////////////////////////////

  /**
   */
  ECELL_ERREST_STEADY_DIFFUSION_T4 (GCELL *gcell);
  
  /**
   */
  void assemble_error (FIELD_SCALAR *field, FIELD_VECTOR *geometry, class PROTO_ERREST<1> *proto);

 private: // class data /////////////////////////////////////////////
  
  static const char TYPE_NAME[];

};

#endif
