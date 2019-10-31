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
#ifndef MAT_MAP_H
# define MAT_MAP_H

#include "mat_mgr.h" 
extern "C" {
#include "tokensP.h"
}

/**
   This class is used to hold instances of materials.
   They may be then retrieved by name.  Protocols would use
   these maps to provide ecells with material handles.
*/
template <typename MAT_CLASS>
class MAT_MAP {

 public:

  MAT_MAP (MAT_MGR *mat_mgr) : _mat_mgr (mat_mgr) {
    _parser = tokens_new_parser ();
  }
  ~MAT_MAP () {
    tokens_delete_parser (_parser);
  }

  /**
     Return the material by name.  It is either found in the map,
     or it is first created and then returned.
  */
  MAT_CLASS *mat (string mattype, string matname) {
    MAT_CLASS *m = dynamic_cast <MAT_CLASS *> (_mat_mgr->mat (mattype, matname));
    if (m == 0) { cerr << "Material " << matname << " of type " << mattype << " not located" << endl; }
    CHECK (m != 0, EXCEPTION_NULL_PTR,;);
    return m;
  }

 private:

  MAT_MGR       *_mat_mgr;
  tokens_parser  _parser;

};

#endif
