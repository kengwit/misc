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
#ifndef BFUN_FE_H
# define BFUN_FE_H

#include <vector>
#include "bfun.h"
#include "enumerator.h"
#include "gcell.h"
#include "ref_fen_src.h"

class GCELL;

class BFUN_FE : public BFUN {

 public: // object functions ////////////////////////////////////////
  
  BFUN_FE (FEN *fen, vector <GCELL *> gcells, BFUN_SET *bfun_set); 
  /**
     Clone a basis function.
  */
  virtual BFUN *clone (BFUN_SET *bfun_set);
  
  // virtual BFUN *clone (FEN*fen,vector<class GCELL*> gcell );
  
  virtual ~BFUN_FE () {}
  virtual  bool refine (REF_CTX *ref_ctx);
  virtual  void unrefine ();

  /** 
   return true if basis function has atleast one inactive parent.
  */ 
  bool has_inactive_parent(BFUN *current_parent);

  vector<GCELL*> collect_child_gcells(FEN *fen);
};

#endif
