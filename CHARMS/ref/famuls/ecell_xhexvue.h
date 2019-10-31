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
#ifndef ECELL_XHEXVUE_H
# define ECELL_XHEXVUE_H

#include "ecell.h"
#include "field.h"
#include "ioP.h"
#include "evalpt.h"
#include "quant_evaluator.h"

/**
Evaluation cell that supports the protocol writing extended HexVue files.
Works together with the protocol PROTO_XHEXVUE.

*/
class ECELL_XHEXVUE : public ECELL {

 public: // declarations ////////////////////////////////////////////

  typedef ECELL_XHEXVUE *(*make_ecell_func) (GCELL *gcell);

 public: // class functions  ////////////////////////////////////////

  /**
     Call this function to register a constructor for the ecell
     based on the gcell type specified, and of given implementation.
  */
  static bool register_make_func (make_ecell_func make, string gcell_type, string implementation);
  static ECELL_XHEXVUE *make_ecell (GCELL *gcell, string implementation);

 public: // object functions ////////////////////////////////////////

  /**
   */
  ECELL_XHEXVUE (GCELL *gcell) : ECELL (gcell) {}
  /**
     Write the HEXVUE record for the ecell to the IO buffer.
     Use the quantity evaluator to compute the values of the desired
     quantities at the nodes or other locations within the cell.
   */
  //  virtual void write (ckit_iobuf_t IOB, QUANT_EVALUATOR &qe) = 0; //ashi
  virtual void write (ckit_iobuf_t IOB, SMART_HANDLE<QUANT_EVALUATOR>qe) = 0;

  // virtual ~ECELL_HEXVUE () {}; //RAW:ashi

  virtual ~ECELL_XHEXVUE () {};
 private: // class objects  /////////////////////////////////////////

  //  static OFACTORY<ECELL_HEXVUE, make_ecell_func, GCELL *> *_ecell_factory; //ashi
 static OFACTORY<ECELL_XHEXVUE, make_ecell_func, GCELL *> *_ecell_factory;

 private: // object data ////////////////////////////////////////////

};



#endif
