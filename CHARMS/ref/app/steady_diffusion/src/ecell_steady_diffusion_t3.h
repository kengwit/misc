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
#ifndef ECELL_STEADY_DIFFUSION_T3_H
# define ECELL_STEADY_DIFFUSION_T3_H

#include "ecell.h"
#include "ecell_steady_diffusion.h"
#include "field.h"
#include "proto_steady_diffusion.h"
#include "trinumintgP.h"
#include "evalpt.h"

class ECELL_STEADY_DIFFUSION_T3 : public ECELL_STEADY_DIFFUSION {

 public: // class functions /////////////////////////////////////////

  static bool register_make_func ();

 public: // object functions ////////////////////////////////////////

  ECELL_STEADY_DIFFUSION_T3 (GCELL *gcell);

  bool assemble_conductivity_matrix (FIELD_SCALAR *phi, LES *les, FIELD_VECTOR *geometry);

  bool assemble_source_terms (FIELD_SCALAR *phi, LES *les, FIELD_VECTOR *geometry);

 private: // class data /////////////////////////////////////////////

  static const char TYPE_NAME[];
  
};

#endif
