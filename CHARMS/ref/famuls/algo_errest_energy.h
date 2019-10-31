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
#ifndef ALGO_ERREST_ENERGY_H
#define ALGO_ERREST_ENERGY_H

#include "fen.h"
#include "algo.h"
#include "gmesh.h"
#include "db.h"
#include "field.h"
#include "proto_errest.h"
#include "algo_errest.h"

template <int NUM_COMPONENTS>
class ALGO_ERREST_ENERGY : public ALGO_ERREST<NUM_COMPONENTS> {
  
 public: // object methods //@@//@@//@@//@@//@@//@@//@@//@@//@@//@@//
  
  ALGO_ERREST_ENERGY (string name, MGR *mgr);
  ~ALGO_ERREST_ENERGY ();
  
  void assemble_error (FIELD_VECTOR *geometry, FIELD<NUM_COMPONENTS> *field);
  void assemble_error (FIELD_VECTOR *geometry, FIELD<NUM_COMPONENTS> *field, double time);
  void accumulate_error (FIELD_VECTOR *geometry, FIELD<NUM_COMPONENTS> *field, double time);
  
  double errval (BFUN *bfun);

 private: // object data //@@//@@//@@//@@//@@//@@//@@//@@//@@//@@//@@

  PROTO_ERREST<NUM_COMPONENTS> *_proto_errest;

};

template <int NUM_COMPONENTS>
ALGO_ERREST_ENERGY<NUM_COMPONENTS>::ALGO_ERREST_ENERGY (string name, MGR *mgr)
  : ALGO_ERREST<NUM_COMPONENTS>  (name, mgr)
{
  string path = "algorithms/errest/" + this->name();
  _proto_errest = 0;
}

template <int NUM_COMPONENTS>
ALGO_ERREST_ENERGY<NUM_COMPONENTS>::~ALGO_ERREST_ENERGY ()
{
  if (_proto_errest) delete _proto_errest;
}

template <int NUM_COMPONENTS>
double 
ALGO_ERREST_ENERGY<NUM_COMPONENTS>::errval (BFUN *bfun)
{
  CHECK (_proto_errest, EXCEPTION_NULL_PTR,;);
  return _proto_errest->errval (bfun);
}

template <int NUM_COMPONENTS>
void
ALGO_ERREST_ENERGY<NUM_COMPONENTS>::assemble_error (FIELD_VECTOR *geometry,
                                                    FIELD<NUM_COMPONENTS> *field) {
  if (_proto_errest) delete _proto_errest;
  _proto_errest = new PROTO_ERREST<NUM_COMPONENTS> (this);
  _proto_errest->assemble_error (geometry, field);
}

template <int NUM_COMPONENTS>
void
ALGO_ERREST_ENERGY<NUM_COMPONENTS>::assemble_error (FIELD_VECTOR *geometry,
                                                    FIELD<NUM_COMPONENTS> *field, double time) {
  assemble_error (geometry, field);
}

template <int NUM_COMPONENTS>
void
ALGO_ERREST_ENERGY<NUM_COMPONENTS>::accumulate_error (FIELD_VECTOR *geometry,
                                                    FIELD<NUM_COMPONENTS> *field, double time) {
  if (!_proto_errest) _proto_errest = new PROTO_ERREST<NUM_COMPONENTS> (this);
  _proto_errest->assemble_error (geometry, field);
}

#endif
