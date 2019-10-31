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
#ifndef CONN_BASE_H
# define CONN_BASE_H

#include "ptr_vector.h"
#include "fen.h"

typedef enum CONN_MANIFOLD_DIM {
  CONN_0_MANIFOLD = 0,
  CONN_1_MANIFOLD = 1,
  CONN_2_MANIFOLD = 2,
  CONN_3_MANIFOLD = 3
}            CONN_MANIFOLD_DIM;

typedef unsigned int CONN_CODE;

class CONN_BASE {

 public: // object pure virtual functions ////////////////////////////

  /**
   */
  CONN_BASE () {}
  virtual ~CONN_BASE () {}

  /**
     Return the number of nodes that this connectivity binds.
   */
  virtual sizet nfens () const = 0;

  /**
     Return the i-th node.
  */
  virtual FEN *fen (sizet i) const = 0;

  /**
     Return the i-th node for modification.
  */
  virtual FEN *& fen (sizet i) = 0;

  /**
     What is the manifold dimension of the conn?
     (Point=0, line=1, surface=2, or solid=3.) 
   */
  virtual CONN_MANIFOLD_DIM manifold_dim () const = 0;

  /**
     Get connectivity code.
  */
  virtual CONN_CODE conn_code () = 0;

  /**
     Return the number of refinement nodes that this connectivity uses
     in a refinement step.
   */
  virtual sizet nreffens () const = 0;
  
  /**
     Get a refinement node.
     The refinement nodes are oriented for the connectivity on
     input.  Therefore, first try to see if there is any topological
     transformation that may be applied to the connectivity on input
     that would lead to a successful match with self.  If that is
     the case, the same topological operations need to be applied
     to the index into vector of refinement nodes indx, and
     ref_fens[transf(indx)] is returned; otherwise, that is if match
     is not possible, null is returned.
  */
  virtual FEN * get_ref_fen (CONN_BASE *conn, std::vector <FEN *> ref_fens, sizet indx) = 0;

  /**
     Clone connectivity.
  */
  virtual CONN_BASE *clone () = 0;

};

#endif
