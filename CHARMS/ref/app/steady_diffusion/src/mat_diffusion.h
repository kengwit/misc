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
#ifndef MAT_DIFFUSION_H
# define MAT_DIFFUSION_H

#include "mat.h"
#include "point.h"
#include "square_fixed_matrix.h"
#include "func.h"

class MAT_DIFFUSION : public MAT {
  
 public: // object methods //////////////////////////////////////////

  MAT_DIFFUSION (DB *db, string name);
  MAT_DIFFUSION (string name);

  /**
     Return the type name.
   */
  virtual string type_name () const  = 0;
  /**
     Return the conductivity coefficient.
     This method is useful in that it makes possible
     to use both isotropic and anisotropic materials
     from within code that can deal with isotropic
     materials only.
  */
  virtual double conductivity (POINT &at) = 0;
  /**
     Return internal heat generation density.
  */
  virtual double internal_heat_generation_density (POINT &at)  = 0;
  /**
     Return value of specified variable.
   */
  double var (string name, POINT &at)  = 0;

};

#endif

