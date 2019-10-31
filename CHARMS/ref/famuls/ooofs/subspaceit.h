/*
                           OOOFS
                 Copyright (C) 2003, Petr Krysl
                 
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/
#ifndef ssit_h
#   define ssit_h

#include <stdio.h>
#include "sparsemtrx.h"
#include "flotarry.h"
#include "logger_stream.h"
#include "logger_stream_mgr.h"

namespace OOOFS {

/** 
Subspace iteration: block inverse power method.
   K y = (omega)^2 M y 
*/
class SubspaceIteration 
{

public :

  SubspaceIteration ();
  ~SubspaceIteration () ;

  /**
  Solve Kx = o^2 Mx. K coefficient matrix, M coefficient matrix, x eigen vector(s), o eigen value(s),
  rtol tolerance, nroot number of required eigenvalues
  */
  bool solve (LOGGER_STREAM &ls,
              SparseMtrx* K, SparseMtrx* M,
              FloatArray* o, FloatMatrix* x, double rtol, int max_iter, int nroot);

private:

  FloatMatrix* ar;
  FloatMatrix* br;
  FloatMatrix*   vec;
  int            n,nc,nsmax,nitem;
  int            solved;

};

} 

#endif









