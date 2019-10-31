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
#ifndef REF_CTX_H
# define REF_CTX_H

#include "ref_fen_src.h"
#include "enumerator.h"
#include <set>

class REF_CTX {

 public: // object methods //////////////////////////////////////////


 public: // object methods //////////////////////////////////////////

  /**
     A refinement context.
     A refinement context holds the incremental changes to the basis
     function set resulting from refinement or derefinement.

     A refinement context holds a set of active
     nodes and nodes to be deactivated.  It takes 
     a field (or perhaps several fields) and refines the geometric
     cells to build up a geometrical structure for basis function refinement.
     The nodes of the basis functions to deactivate are also listed
     in the refinement context as inactive nodes.  
  */
  REF_CTX (GMESH *gmesh) {
    _gmesh = gmesh;
    _ref_fen_src = new REF_FEN_SRC ();
    _activated_fens.clear ();
    _deactivated_fens.clear ();
    _true_hierarchical = true;
  }
  ~REF_CTX () {
    if (_ref_fen_src) delete _ref_fen_src;
  }
  /**
   */
  bool true_hierarchical () { return _true_hierarchical; }
  /**
   */
  void set_true_hierarchical (bool true_hierarchical) { _true_hierarchical = true_hierarchical; }

  /**
     Get a refinement node --- see REF_FEN_SRC.
  */
  FEN *get_ref_fen (CONN_BASE *conn, sizet indx, vector <POINT> ref_locs) {
    FEN_MAKER fen_maker (_gmesh, ref_locs);
    return _ref_fen_src->get_ref_fen (conn, indx, &fen_maker);
  }
 
 private: // object data ////////////////////////////////////////////
  GMESH                *_gmesh;
  REF_FEN_SRC          *_ref_fen_src;
  set <FEN *>           _activated_fens;
  set <FEN *>           _deactivated_fens;
  bool                  _true_hierarchical;
  
 private: // object methods ////////////////////////////////////////
  
};

#endif
