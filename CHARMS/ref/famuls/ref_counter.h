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
#ifndef REF_COUNTER_H
#  define REF_COUNTER_H

#include "famexception.h"

/**
   Reference counting object.
*/
class REF_COUNTER {
  
 public: // object methods ///////////////////////////////////////////
  
  /**
     Constructor.
  */
  inline REF_COUNTER (void) :	_count (1) {}
  
  /**
     Add reference to the counter.  The counter returns
     false if the counter is no longer valid.
  */
  inline bool add_ref (void) {
    if (_count > 0) {
      ++_count;
      return true;
    } else
      return false;
  }
  
  /**
     Remove one reference.  The object that pointed at this
     counter should not use it anymore.
     Returns:
     o false (nothing needs to be done);
     o true: delete both the counter
     and the object whose uses are being counted.
  */
  inline bool remove_ref (void) {
    CHECK (_count > 0, EXCEPTION_BAD_ACCESS,;);
    --_count;
    if (_count == 0) {
      return true;
    }  else
      return false;
  }
  
 private: // private methods //////////////////////////////////////////

  /**
     Calls prevented to avoid copying.
  */
  REF_COUNTER (const REF_COUNTER&);
  REF_COUNTER& operator= (const REF_COUNTER&);
  
 private: // object data /////////////////////////////////////////////
  
  unsigned long	_count;
  
};

#endif
