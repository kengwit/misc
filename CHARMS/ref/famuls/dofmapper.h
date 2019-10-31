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
#ifndef DOFMAPPER_H
# define DOFMAPPER_H

#include "dofmapper_base.h"

class FEN;

template <int NUM_COMPONENTS>
class DOFMAPPER : public DOFMAPPER_BASE {

 public: // object functions ////////////////////////////////////////

  DOFMAPPER (class FEN *fen,
             BFUN_DOFPARAM_PAIR_ID dofparam_id, bool constrained[NUM_COMPONENTS])
    : DOFMAPPER_BASE::DOFMAPPER_BASE (fen, dofparam_id) {
    for (int j = 0; j < NUM_COMPONENTS; j++) {
      _eqnums[j] = DOFMAPPER_BASE::invalid_eqnum;
      _constrained[j] = constrained[j];
    }
  }
  /**
     Set equation numbers. The dofmapper object knows which equations are constrained
     so remember, constrained equations will be set to DOFMAPPER_BASE::invalid_eqnum.
  */
  void set_eqnums (int eqnums[NUM_COMPONENTS]) {
    for (int j = 0; j < NUM_COMPONENTS; j++) {
      if (_constrained[j]) _eqnums[j] = DOFMAPPER_BASE::invalid_eqnum;
      else                 _eqnums[j] = eqnums[j];
    }
  }
  /**
     How many equations?
  */
  int neqns () const { return NUM_COMPONENTS; }
  /**
     Get equation numbers.  If an equation is constrained,
     it will be returned as invalid_eqnum.
  */
  void eqnums (int eqnums[NUM_COMPONENTS]) const {
    for (int j = 0; j < NUM_COMPONENTS; j++) {
      if (_constrained[j]) eqnums[j] = DOFMAPPER_BASE::invalid_eqnum;
      else                 eqnums[j] = _eqnums[j];
    }
  }
  /**
     Get constraint flags.  If an equation is constrained,
     the flag will be true; otherwise false.
  */
  void constrained (bool c[NUM_COMPONENTS]) const {
    for (int j = 0; j < NUM_COMPONENTS; j++) {c[j] = _constrained[j];}
  }
  
 private: // object data ////////////////////////////////////////////

  int            _eqnums[NUM_COMPONENTS];
  bool           _constrained[NUM_COMPONENTS];

};

#endif
