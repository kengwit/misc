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
#ifndef ECELL_ERREST_H
# define ECELL_ERREST_H

#include "ecell.h"
#include "field.h"
#include "ofactory.h"

template <int NUM_COMPONENTS> class PROTO_ERREST;

template <int NUM_COMPONENTS>
class ECELL_ERREST : public ECELL {

 public: // declarations ////////////////////////////////////////////

  typedef ECELL_ERREST *(*make_ecell_func) (GCELL *gcell);

 public: // class functions  ////////////////////////////////////////

  /**
     Call this function to register a constructor for the ecell
     based on the gcell type specified, and of given implementation.
  */
  static bool register_make_func (make_ecell_func make, string gcell_type, string implementation);
  static ECELL_ERREST *make_ecell (GCELL *gcell, string implementation);

 public: // object functions ////////////////////////////////////////

  ECELL_ERREST () : ECELL(0) {}
  /**
   */
  ECELL_ERREST (GCELL *gcell) : ECELL (gcell) {}
  /**
     Attach the ecell to the protocol.  The protocol lets the ecell
     retrieve any data needed by calling this function.
   */
  virtual void attach (class PROTO_ERREST<NUM_COMPONENTS> *proto, string ggpath) {}
  /**
   */
  virtual void assemble_error (FIELD<NUM_COMPONENTS> *field,
                               FIELD_VECTOR *geometry,
                               class PROTO_ERREST<NUM_COMPONENTS> *proto) = 0;
  /**
   */
  virtual ~ECELL_ERREST () {};

 private: // class objects  /////////////////////////////////////////

  static OFACTORY<ECELL_ERREST, make_ecell_func, GCELL *> *_ecell_factory;



};


template <int NUM_COMPONENTS>
OFACTORY<ECELL_ERREST<NUM_COMPONENTS>,
         typename ECELL_ERREST<NUM_COMPONENTS>::make_ecell_func,
         GCELL*> *ECELL_ERREST<NUM_COMPONENTS>::_ecell_factory = 0;

template <int NUM_COMPONENTS>
bool
ECELL_ERREST<NUM_COMPONENTS>::register_make_func (ECELL_ERREST<NUM_COMPONENTS>::make_ecell_func make,
                                                  string gcell_type, string implementation)
{
  if (_ecell_factory == 0) {
    ECELL_ERREST<NUM_COMPONENTS>::_ecell_factory
      = new OFACTORY<ECELL_ERREST<NUM_COMPONENTS>,
      typename ECELL_ERREST<NUM_COMPONENTS>::make_ecell_func, GCELL* >;
  }
  return (ECELL_ERREST<NUM_COMPONENTS>::_ecell_factory)->register_make_func (make, gcell_type, implementation);
}


template <int NUM_COMPONENTS>
ECELL_ERREST<NUM_COMPONENTS> *
ECELL_ERREST<NUM_COMPONENTS>::make_ecell (GCELL *gcell, string implementation)
{
  CHECK (ECELL_ERREST<NUM_COMPONENTS>::_ecell_factory, EXCEPTION_NULL_PTR,;);
  return (ECELL_ERREST<NUM_COMPONENTS>::_ecell_factory)->make (gcell, gcell->type_name (), implementation);
}

#endif
