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
#ifndef FEN_MAKER_H
# define FEN_MAKER_H

#include "fen.h"
#include "gmesh.h"

class FEN_MAKER {

 public:

  FEN_MAKER (GMESH *gmesh, vector <POINT> ref_locs) {
    _gmesh = gmesh; _ref_locs = ref_locs;
  }    
  
  FEN *operator () (int i) {
    FEN *fen = new FEN (_gmesh->max_fen_id()+1, _ref_locs[i], _gmesh);
    return fen;
  }
  
 private:

  GMESH         *_gmesh;
  vector <POINT> _ref_locs;
  
};

#endif
