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
#include "fixed_vector.h"

FIXED_VECTOR<3> cross_prod( const FIXED_VECTOR<3>& u, const FIXED_VECTOR<3>& v ) {
  FIXED_VECTOR<3> result;
  result(0) =   u[1]*v(2) - u[2]*v(1);
  result(1) =   u[2]*v(0) - u[0]*v(2);
  result(2) =   u[0]*v(1) - u[1]*v(0);
  return result;
}

#if 0
extern "C" {
#include "matrix.h"
#include "randP.h"
}
#include "famuls.h"
int
main ()
{
  const double eps = 1e-7;
  const int dim = 18;
#define MAKE_VECTS(name)                                                \
  VEC *name = v_get (dim); FIXED_VECTOR<dim> v##name;                   \
  { for (int j = 0; j < dim; j++) name->ve[j] = v##name[j] = rand_0_to_1 (); }
#define COMP_VECTS(name)                                                                \
  {                                                                                     \
    bool result = true;                                                                 \
    for (int j = 0; j < dim; j++)   {                                                   \
      if (fabs (name->ve[j] - v##name[j]) > eps) {                                     \
        cerr << "mismatch: " << STRINGIZE_(name) << "[" << j << "] =" << name->ve[j] << \
          " vs. " << STRINGIZE_(v##name) << "[" << j << "] =" << v##name[j] << endl;    \
      result = false;                                                                   \
      }                                                                                 \
    }                                                                                   \
    if (result) cerr << "Correct" << endl;                                              \
    else        cerr << "Error" << endl;                                                \
  }

  MAKE_VECTS (a);
  MAKE_VECTS (b);
  MAKE_VECTS (c);
  MAKE_VECTS (d);
  MAKE_VECTS (e);
  MAKE_VECTS (f);
  MAKE_VECTS (g);

  vf = vg; f = v_copy (g, f);  COMP_VECTS (f);

  vf += vb; f = v_add (f, b, f); COMP_VECTS (f);

  if (fabs (v_norm1 (a) - va.l1_norm ()) > eps) cerr << "Error" << endl;
  else                                          cerr << "Correct" << endl;

  if (fabs (v_norm2 (a) - va.l2_norm ()) > eps) cerr << "Error" << endl;
  else                                          cerr << "Correct" << endl;
  COMP_VECTS (a);

  if (fabs (v_norm2 (b) - vb.l2_norm ()) > eps) cerr << "Error" << endl;
  else                                          cerr << "Correct" << endl;

  if (fabs (v_norm_inf (b) - vb.linfinity_norm ()) > eps) cerr << "Error" << endl;
  else                                                    cerr << "Correct" << endl;

  va[dim-1] = a->ve[dim-1]; COMP_VECTS (a);
  b->ve[dim-1] = vb(dim-1); COMP_VECTS (b);
  vb += vc; v_add (b, c, b); COMP_VECTS (b);
  vb -= vd; v_sub (b, d, b); COMP_VECTS (b);
  vb = 3.0; v_ones (b); sv_mlt (3.0, b, b); COMP_VECTS (b);
  vb.zero (); v_zero (b); COMP_VECTS (b);

  if (fabs (ve*vg - in_prod(e,g)) > eps) cerr << "Error" << endl;
  else                                   cerr << "Correct" << endl;

  if (fabs (vc*vg - in_prod(c,g)) > eps) cerr << "Error" << endl;
  else                                   cerr << "Correct" << endl;

  ve.add (-5, vf); v_mltadd (e, f, -5, e); COMP_VECTS (e);

  ve.add (-5, vf); v_mltadd (e, f, -5, e); COMP_VECTS (e);

  va.scale (-0.533); sv_mlt (-0.533, a, a); COMP_VECTS (a);

  va.assign (-1, vc); v_copy (c, a); sv_mlt (-1, a, a); COMP_VECTS (a);

  vd.assign (0.22, va, -0.3, ve); v_zero (d); v_mltadd (d, a, 0.22, d); v_mltadd (d, e, -0.3, d); COMP_VECTS (d);

  vd = va - vb; v_sub (a, b, d); COMP_VECTS (d);

  vd = va + vd - ve + vf.scale (vd * vc);
  VEC *tmp = v_get (dim);
  for (int j = 0; j < dim; j++) {
    tmp->ve[j] = a->ve[j] + d->ve[j] - e->ve[j] + in_prod (d, c) * f->ve[j];
  }
  d = v_copy (tmp, d);
  COMP_VECTS (d);
}
#endif
