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
#ifndef ALGO_ERREST_GRADATION_H
#define ALGO_ERREST_GRADATION_H

#include "fen.h"
#include "algo.h"
#include "gmesh.h"
#include "db.h"
#include "field.h"
#include "algo_errest.h"
#include "proto_errest.h"
#include "algo_gradation.h"

template <int NUM_COMPONENTS>
class ALGO_ERREST_GRADATION : public ALGO_ERREST<NUM_COMPONENTS> {
  
 public: // object methods //@@//@@//@@//@@//@@//@@//@@//@@//@@//@@//
  
  ALGO_ERREST_GRADATION (string name, MGR *mgr, GMESH *gmesh);
  ALGO_ERREST_GRADATION (string name, MGR *mgr, GMESH *gmesh, double t);
  ~ALGO_ERREST_GRADATION ();

  void assemble_error (FIELD_VECTOR *geometry, FIELD<NUM_COMPONENTS> *field);
  void assemble_error (FIELD_VECTOR *geometry, FIELD<NUM_COMPONENTS> *field,double ttime);

  double errval (BFUN *bfun);

  void write_mesh_size_to_name (string filename);

 private: // object data //@@//@@//@@//@@//@@//@@//@@//@@//@@//@@//@@

  map <BFUN *, double> _bfun_error_map;
  ALGO_GRADATION      *_algo_gradation;

 private:
  
  map <BFUN *, double> do_assemble (FIELD<NUM_COMPONENTS> *field, FIELD_VECTOR *geometry, ALGO_GRADATION *_algo_gradation);

};

template <int NUM_COMPONENTS>
ALGO_ERREST_GRADATION<NUM_COMPONENTS>::ALGO_ERREST_GRADATION (string name, MGR *mgr, GMESH *gmesh)
  : ALGO_ERREST<NUM_COMPONENTS>  (name, mgr)
{
  double t = 0.0;
  string path = "algorithms/errest/" + this->name();
  _algo_gradation = new ALGO_GRADATION (this->name(), mgr, gmesh,t);
}

template <int NUM_COMPONENTS>
ALGO_ERREST_GRADATION<NUM_COMPONENTS>::ALGO_ERREST_GRADATION (string name, MGR *mgr, GMESH *gmesh,double t)
  : ALGO_ERREST<NUM_COMPONENTS>  (name, mgr)
{
  string path = "algorithms/errest/" + this->name();
  _algo_gradation = new ALGO_GRADATION (this->name(), mgr, gmesh,t);
}

template <int NUM_COMPONENTS>
ALGO_ERREST_GRADATION<NUM_COMPONENTS>::~ALGO_ERREST_GRADATION ()
{
  if (_algo_gradation) delete _algo_gradation;
}

template <int NUM_COMPONENTS>
double 
ALGO_ERREST_GRADATION<NUM_COMPONENTS>::errval (BFUN *bfun)
{
  map <BFUN *, double>::const_iterator i = _bfun_error_map.find (bfun);
  if (i != _bfun_error_map.end()) return i->second;
  else                            return 0;
}

template <int NUM_COMPONENTS>
map <BFUN *, double>
ALGO_ERREST_GRADATION<NUM_COMPONENTS>::do_assemble (FIELD<NUM_COMPONENTS> *field,
                                                    FIELD_VECTOR *geometry,
                                                    ALGO_GRADATION *_algo_gradation)
{
  map <BFUN *, double> m;
  for (sizet j = 0; j < field->npairs (); j++) {
    FIELD_PAIR<NUM_COMPONENTS> *fp = field->field_pair (j);
    BFUN *bfun = fp->bfun();
    POINT _p = bfun->fen()->ref_loc();
    G3D_box_t box;
    G3D_loc_t p;
    p.x = _p(0); p.y = _p(1); p.z = _p(2);
    G3D_INIT_BOX (box, p);
    for (sizet k = 0; k < bfun->ngcells(); k++) {
      GCELL *gcell = bfun->gcell(k);
      for (sizet l = 0; l < gcell->conn()->nfens(); l++) {
        _p = gcell->conn()->fen(l)->ref_loc();
        p.x = _p(0); p.y = _p(1); p.z = _p(2);
        G3D_UPDT_BOX (box, p);
      }
    }
    double box_size = (G3D_BOX_RANGE (box, x) +
                       G3D_BOX_RANGE (box, y) +
                       G3D_BOX_RANGE (box, z)) / 3 / 2;
    EVALPT evalpt (field, bfun->gcell(0), bfun->fen());
    double mesh_size = _algo_gradation->mesh_size (evalpt);
    if (mesh_size <= 0) {
      cerr << "Mesh size " << mesh_size << " at " << bfun->fen()->ref_loc() << endl;
      mesh_size = 0.0000001;
    }
    double error = box_size / mesh_size;
    m.insert (map <BFUN *, double>::value_type (bfun, error));
  }
  return m;
}

template <int NUM_COMPONENTS>
void
ALGO_ERREST_GRADATION<NUM_COMPONENTS>::assemble_error (FIELD_VECTOR *geometry,
                                                       FIELD<NUM_COMPONENTS> *field) {
   
   _bfun_error_map = do_assemble (field, geometry, _algo_gradation);
   
}

template <int NUM_COMPONENTS>
void
ALGO_ERREST_GRADATION<NUM_COMPONENTS>::assemble_error (FIELD_VECTOR *geometry,
                                                       FIELD<NUM_COMPONENTS> *field, double ttime) {

  _algo_gradation->init_mesh_size(ttime);
  _bfun_error_map = do_assemble (field, geometry, _algo_gradation);
}

template <int NUM_COMPONENTS>
void
ALGO_ERREST_GRADATION<NUM_COMPONENTS>::write_mesh_size_to_name (string filename)
{
  _algo_gradation->write_mesh_size_to_name (filename);
}

#endif
