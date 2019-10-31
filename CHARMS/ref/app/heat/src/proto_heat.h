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
#ifndef PROTO_HEAT_H
# define PROTO_HEAT_H

/*
  Transient heat equation: M(T) T^dot + K(T) T = F(T)

  rho C T^dot = d/dx_i (k_ij dT/dx_j) + Q

  C = specific heat, k_ij = cartesian components of the conductivity
  tensor (symmetric), Q = internal heat generation
  
  with bc's

  T = f_T (x, t) on Gamma_T

  -(k_ij dT/dx_j) n_i = q_a + q_c + q_r = f_q (x, t) on Gamma_q

  q_a = applied flux
  q_c = convective flux, q_c = h_c(x, T, t) (T - T_c)
  q_r = radiative flux, q_r = h_r(x, T, t) (T - T_r)

 */

#include <list>
#include "ecell_heat.h"

class PROTO_HEAT {

 public:

  PROTO_HEAT ();
  bool PROTO_HEAT::assemble_conductivity_matrix ();
  bool PROTO_HEAT::assemble_capacitance_matrix ();
  bool PROTO_HEAT::assemble_source_terms ();
  
 private:

  std::list <ECELL_HEAT *> _ecells;
  
};

#endif
