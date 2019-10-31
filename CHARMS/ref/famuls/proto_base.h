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
#ifndef PROTO_BASE_H
# define PROTO_BASE_H

#include <list>
#include "algo.h"
#include "gmesh.h"
#include "ecell.h"
#include "field.h"
#include "quant_evaluator.h"
#include "enumerator.h"
extern "C" {
#include "ioP.h"
}

/**
This is the base class for the protocol.  
It is in the form of an interface. 
*/
class PROTO_BASE {

 public: //	class declarations //////////////////////////////////////

	 /**
	 This class needs to be overridden 
	 a protocol that inherits from PROTO_BASE. 	 
	 */
  class ECELL_ENUMERATOR  {
    public :
    ECELL_ENUMERATOR (PROTO_BASE *proto) {};
    virtual ECELL *next() = 0;
    virtual void reset()  = 0;
   };

 public: //	object functions ////////////////////////////////////////

  PROTO_BASE (){};
  virtual ~PROTO_BASE () {}

  /**
  Enumerate evaluation cells.  
  */
 
  virtual ECELL_ENUMERATOR *ecell_enumerator () = 0;
  /**
  Make a quantity evaluator for the given ecell.
  */
  
  virtual SMART_HANDLE<QUANT_EVALUATOR> quant_evaluator (ECELL *e,list<string> quants) = 0;
  
};

#endif
