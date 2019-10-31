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
#include "famuls.h"
#include "ecell_wave_traction_q4.h"
#include "evalpt.h"
#include "square_fixed_matrix.h"
#include "gcell_surf_q4.h"
#include "proto_wave.h"

const char ECELL_WAVE_TRACTION_Q4::_type_name[] = "wave_traction_q4";

ECELL_WAVE_TRACTION_Q4::ECELL_WAVE_TRACTION_Q4 (GCELL *gcell)
  : ECELL_WAVE::ECELL_WAVE (gcell)
{
}



bool 
ECELL_WAVE_TRACTION_Q4::assem_tractions_load ()
{
  FIELD_VECTOR *geometry = _proto->geometry();
  FIELD_VECTOR *u = _proto->u();
  LOAD_ON_GCELL_TRACTION *load = this->traction_load();
  if (load) {
    GAUSS_2d_integration_rule_t *gr = gauss_rule_2d (2);
    for (int qp = 0; qp < gr->npoints; qp++) {
      POINT param_loc (gr->ip[qp].xi,gr->ip[qp].eta , 0);
      EVALPT evalpt (u, gcell (), param_loc);
      evalpt.eval ();
      int nbfuns = evalpt.nbfuns () ;
      POINT at = geometry->evaluate (geometry, this->gcell(), param_loc);
      POINT4 at4; at4(0) = at(0);at4(1) = at(1);at4(2) = at(2);at4(3) = _proto->time();
      FUNC<3> &f = load->func();
      FIXED_VECTOR<3> traction = f(at4);        
      //FIXED_VECTOR<3> traction = load->traction(at4);
      FIXED_VECTOR<3> tangent_vector(0);
      FIXED_VECTOR<3> tangent_vector1(0);
      for (int I = 0; I < nbfuns; I++) {
        double Nx_I = evalpt.N_x(I);
        double Ny_I = evalpt.N_y(I);
        FEN *fen =  evalpt.fen(I);
        BFUN_DOFPARAM_PAIR_ID dpid = geometry->dofparam_id (fen);
        FIELD_PAIR<3> *fp = geometry->field_pair (dpid);
        FIXED_VECTOR<3> v = fp->dofparam ();
        tangent_vector.add (Nx_I, v);
        tangent_vector1.add (Ny_I, v);
      }
      FIXED_VECTOR<3> result = cross_prod(tangent_vector,tangent_vector1);
      double fact = gr->ip[qp].Wgt * result.l2_norm();
    
      // now loop over all nodes 
      for (int I = 0; I < nbfuns; I++) {
        double N_I = evalpt.N (I);
        double values[3];
        values[0] = traction(0)*fact*N_I;
        values[1] = traction(1)*fact*N_I;
        values[2] = traction(2)*fact*N_I;
        _proto->assemble_f(evalpt.dofparam_id(I),values);
      }
    }
    gauss_free_2d (gr);
  }
  return true;
}

static ECELL_WAVE *
make (GCELL *gcell)
{
  return (new ECELL_WAVE_TRACTION_Q4 (gcell));
}

bool
ECELL_WAVE_TRACTION_Q4::register_make_func ()
{
  return ECELL_WAVE::register_make_func (make, string (GCELL_SURF_Q4::TYPE_NAME),
                                               string (_type_name));
}

