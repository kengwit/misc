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
#ifndef EVALBUF_H
# define EVALBUF_H

#include "fen.h"
#include "bfun.h"
#include "bfun_dofparam_pair_id.h"
#include "field_base.h"
#include "general_fixed_vector.h"

class EVALBUF {

 private: // object internal definitions ////////////////////////////

  class ROW {
  public:
    ROW () {
      dofparam_id = INVALID_BFUN_DOFPARAM_PAIR_ID;
      fen = 0;
      bfun = 0;
      N = 0;
      for (int j = 0; j < 3; j++) N_der[j] = 0;
    }
    BFUN_DOFPARAM_PAIR_ID       dofparam_id;
    FEN                        *fen;
    BFUN                       *bfun;
    double                      N;
    double                      N_der[3]; // largest possible
  };
  
 public: // object methods //////////////////////////////////////////

  EVALBUF (FIELD_BASE *field, FIELD_BASE::PAIRS_ITERATORS pi) : _rows () {
    CHECK(pi.first != pi.second, EXCEPTION_BAD_ACCESS,;); // no active pairs?
    int n = 0;
    for (set <BFUN_DOFPARAM_PAIR_ID>::const_iterator i = pi.first; i != pi.second; i++) { n++; }
    _rows.size(n);
    _rows.assign(ROW());
    int j = 0;
    for (set <BFUN_DOFPARAM_PAIR_ID>::const_iterator i = pi.first; i != pi.second; i++) { 
      _rows(j).dofparam_id = *i;
      _rows(j).bfun = field->bfun (_rows(j).dofparam_id);
      _rows(j).fen  = _rows(j).bfun->fen();
      j++;
    }
  }
  ~EVALBUF () {}
  ROW &operator () (const sizet j) {
    return _rows(j);
  }
  sizet size () { return _rows.size (); }

 private: // object data ////////////////////////////////////////////

  GENERAL_FIXED_VECTOR <ROW> _rows;

};

#endif
