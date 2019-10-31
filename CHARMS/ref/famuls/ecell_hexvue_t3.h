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
#ifndef ECELL_HEXVUE_T3_H
# define ECELL_HEXVUE_T3_H

#include "ecell_hexvue.h"
#include "gcell_surf_t3.h"

template <int NUM_COMPONENTS>
class ECELL_HEXVUE_T3 : public ECELL_HEXVUE<NUM_COMPONENTS> {

 public: // declarations ////////////////////////////////////////////

  static bool register_make_func ();

 public: // class functions  ////////////////////////////////////////

  ECELL_HEXVUE_T3 (GCELL *gcell)
    : ECELL_HEXVUE<NUM_COMPONENTS>::ECELL_HEXVUE (gcell) {}

  /**
     Write the HEXVUE record for the ecell to the IO buffer.
     The cell will have the geometry given on input, and the field
     to display on the geometry is also specified.  To evaluate the
     field, the geometry eval_on_geometry is used.
  */
  void write (ckit_iobuf_t IOB,
              FIELD_VECTOR *display_on_geometry,
              FIELD_VECTOR *eval_on_geometry,
              FIELD<NUM_COMPONENTS> *field)
    {
      double L[3][3] = {{0,0,1}, {1,0,0}, {0,1,0}};
      const int NFENS = 3;
      FIXED_VECTOR<3> x[NFENS];
      FIXED_VECTOR<NUM_COMPONENTS> v[NFENS];
      
      for (int pt = 0; pt < NFENS; pt++) {
        // evaluate the geometry
        POINT param_loc (L[pt][0], L[pt][1], 0);
        EVALPT gevalpt (display_on_geometry, this->gcell(), param_loc);
        gevalpt.eval ();
        {
          const int nbfuns = gevalpt.nbfuns ();
          x[pt] = 0;
          for (int i = 0; i < nbfuns; i++) {
            double N_i = gevalpt.N (i);
            BFUN *bfun = gevalpt.bfun (i);
            FEN *fen = bfun->fen ();
            BFUN_DOFPARAM_PAIR_ID dpid = display_on_geometry->dofparam_id (fen);
            FIELD_PAIR<3> *fp = display_on_geometry->field_pair (dpid);
            FIXED_VECTOR<3> dofparam = fp->dofparam ();
            x[pt].add(N_i, dofparam);
          }
        }
        // evaluate the field
        EVALPT fevalpt (field, this->gcell(), param_loc);
        fevalpt.eval ();
        {
          const int nbfuns = fevalpt.nbfuns ();
          v[pt] = 0;
          for (int i = 0; i < nbfuns; i++) {
            double N_i = fevalpt.N (i);
            BFUN *bfun = fevalpt.bfun (i);
            FEN *fen = bfun->fen ();
            BFUN_DOFPARAM_PAIR_ID dpid = field->dofparam_id (fen);
            FIELD_PAIR<NUM_COMPONENTS> *fp = field->field_pair (dpid);
            FIXED_VECTOR<NUM_COMPONENTS> dofparam = fp->dofparam ();
            v[pt].add(N_i, dofparam);
          }
        }
      } // for all pt
      write_hexvue_tri3_ascii (IOB, x, v);
    }
  
 private: // object data ////////////////////////////////////////////
  
 private: // class members //////////////////////////////////////////
  
  
  int write_hexvue_tri3_ascii (ckit_iobuf_t IOB,
                               POINT p[3],
                               FIXED_VECTOR<NUM_COMPONENTS> values[3]
                               )
    {
      int j, i;
      char buf[128];
      
#define RUN_IOP(action)                                         \
   {                                                            \
     if (!action) {                                             \
        fprintf (stderr, "Error: %s\n", ((char *)ckit_io_error_explanation ()));     \
     }                                                          \
   }
      RUN_IOP (ckit_io_put_line (IOB, "T3"));
      for (j = 0; j < 3; j++) {
        sprintf (buf, "%g %g %g", p[j](0), p[j](1), p[j](2));
        RUN_IOP (ckit_io_put_line (IOB, buf));
        for (i = 0; i < NUM_COMPONENTS; i++) {
          sprintf (buf, "%g", values[j](i));
          RUN_IOP (ckit_io_put_line (IOB, buf));
        }      
      }
      
      return 1;
    }
  
};

template <int NUM_COMPONENTS>
static ECELL_HEXVUE<NUM_COMPONENTS> *
make (GCELL *gcell)
{
  return (new ECELL_HEXVUE_T3<NUM_COMPONENTS> (gcell));
}

template <int NUM_COMPONENTS>
bool
ECELL_HEXVUE_T3<NUM_COMPONENTS>::register_make_func ()
{
  return ECELL_HEXVUE<NUM_COMPONENTS>::register_make_func (make,
                                                           string (GCELL_SURF_T3::TYPE_NAME),
                                                           string ("default"));
}

#endif
