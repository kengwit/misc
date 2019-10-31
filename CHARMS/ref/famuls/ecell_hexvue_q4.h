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
#ifndef ECELL_HEXVUE_Q4_H
# define ECELL_HEXVUE_Q4_H

#include "ecell_hexvue.h"
#include "gcell_surf_q4.h"

template <int NUM_COMPONENTS>
class ECELL_HEXVUE_Q4 : public ECELL_HEXVUE<NUM_COMPONENTS> {

 public: // declarations ////////////////////////////////////////////

  static bool register_make_func ();

 public: // class functions  ////////////////////////////////////////

  ECELL_HEXVUE_Q4 (GCELL *gcell)
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
      double xi[4][2] = {{-1,-1}, {+1,-1}, {+1,+1}, {-1,+1}};
      const int NFENS = 4;
      FIXED_VECTOR<3> x[NFENS];
      FIXED_VECTOR<NUM_COMPONENTS> v[NFENS];
      
      for (int pt = 0; pt < NFENS; pt++) {
        POINT param_loc (xi[pt][0], xi[pt][1], 0);
        // evaluate the geometry
        x[pt] = display_on_geometry->evaluate (eval_on_geometry, this->gcell(), param_loc);
        // evaluate the field
        v[pt] =  field->evaluate (eval_on_geometry, this->gcell(), param_loc);
      } // for all pt
      write_hexvue_quad4_ascii (IOB, x, v);
    }
  
  
 private: // object data ////////////////////////////////////////////

 private: // class members //////////////////////////////////////////
  
  int write_hexvue_quad4_ascii (ckit_iobuf_t IOB,
                                POINT p[4],
                                FIXED_VECTOR<NUM_COMPONENTS> values[4])
    {
      int j, i;
      char buf[128];
      
#undef RUN_IOP
#define RUN_IOP(action)                                         \
   {                                                            \
     if (!action) {                                             \
        fprintf (stderr, "Error: %s\n", ((char *)ckit_io_error_explanation ()));     \
     }                                                          \
   }
      RUN_IOP (ckit_io_put_line (IOB, "Q4"));
      for (j = 0; j < 4; j++) {
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
  return (new ECELL_HEXVUE_Q4<NUM_COMPONENTS> (gcell));
}

template <int NUM_COMPONENTS>
bool
ECELL_HEXVUE_Q4<NUM_COMPONENTS>::register_make_func ()
{
  return ECELL_HEXVUE<NUM_COMPONENTS>::register_make_func (make,
                                                           string (GCELL_SURF_Q4::TYPE_NAME),
                                                           string ("default"));
}

#endif
