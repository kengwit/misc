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
#ifndef ALGO_XMGR_H
#define ALGO_XMGR_H

#include "fen.h"
#include "algo.h"
#include "gmesh.h"
#include "db.h"
#include "field.h"

class ALGO_XMGR : public ALGO {
  
 public: // object methods //////////////////////////////////////////
  
  ALGO_XMGR (string name, MGR *mgr);
  ~ALGO_XMGR () {}

  void run (string label, FIELD_SCALAR *field, FIELD_VECTOR *geometry);
  
 private: // object data ////////////////////////////////////////////
  
  sizet _set;

  class PAIR_COMPARATOR {
  public:
    bool operator() (FIELD_PAIR<1> *p0, FIELD_PAIR<1> *p1) const {
      return (p0->bfun()->fen()->ref_loc()(0) < p1->bfun()->fen()->ref_loc()(0));
    }
  };

};

#endif
