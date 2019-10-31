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
#ifndef CONN_REF_H
# define CONN_REF_H

#include <vector>
#include "conn.h"

class CONN_REF {

 public: // object functions ////////////////////////////////////////

  /**
   */
  CONN_REF (CONN_BASE *conn) { _conn = conn->clone (); }
  ~CONN_REF () { delete _conn; }
  /**
     Get refinement node.  The task is delegated
     to the underlying connectivity.
  */
  FEN *get_ref_fen (CONN_BASE *conn, sizet indx);

  /**
     Refine the conn_ref by generating refinement nodes.
  */
  void refine (class FEN_MAKER *fen_maker);
  
 private: // object data ////////////////////////////////////////////

  CONN_BASE                 *_conn;
  std::vector <FEN *>        _ref_fens;
  
};

#endif
