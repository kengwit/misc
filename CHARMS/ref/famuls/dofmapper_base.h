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
#ifndef DOFMAPPER_BASE_H
# define DOFMAPPER_BASE_H

#include "bfun_dofparam_pair_id.h"

class FEN;

class DOFMAPPER_BASE {

 public: // class members ///////////////////////////////////////////
  
  static const int invalid_eqnum = -1;

 public: // object functions ////////////////////////////////////////

  DOFMAPPER_BASE (class FEN *fen, BFUN_DOFPARAM_PAIR_ID dofparam_id) {
    _fen = fen;
    _dofparam_id = dofparam_id;
  }
  /**
     Get the fen.
  */
  class FEN *fen () const { return _fen; }
  /**
     Get the dofparam id.
  */
  BFUN_DOFPARAM_PAIR_ID dofparam_id () const { return _dofparam_id; }

 public: // virtual object functions ////////////////////////////////
  
  virtual void set_eqnums (int eqnums[]) = 0;
  virtual int neqns () const = 0;
  virtual void eqnums (int eqnums[]) const = 0;
  virtual void constrained (bool c[]) const = 0;
  
 private: // object data ////////////////////////////////////////////

  class FEN *_fen;
  BFUN_DOFPARAM_PAIR_ID _dofparam_id;

};

#endif
