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
#ifndef ECELL_HEXVUE_T10_H
# define ECELL_HEXVUE_T10_H

#include "ecell_hexvue.h"
#include "gcell_solid_t10.h"

template <int NUM_COMPONENTS>
class ECELL_HEXVUE_T10 : public ECELL_HEXVUE<NUM_COMPONENTS> {

 public: // declarations ////////////////////////////////////////////

  static bool register_make_func ();

 public: // class functions  ////////////////////////////////////////

  ECELL_HEXVUE_T10 (GCELL *gcell)
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
      const int NFENS = this->gcell()->conn()->nfens();
      FIXED_VECTOR<3> x[15];
      FIXED_VECTOR<NUM_COMPONENTS> v[15];
      
      for (int pt = 0; pt < 15; pt++) {
        // evaluate the geometry
        POINT param_loc (g3dtet10_vertex_param_coords[pt][0],
                         g3dtet10_vertex_param_coords[pt][1],
                         g3dtet10_vertex_param_coords[pt][2]);
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
      const int hexi[4][8] = {
        {10,5,2,6,14,12,9,13},
        {4,1,5,10,11,8,12,14},
        {0,4,10,6,7,11,14,13},
        {7,11,14,13,3,8,12,9}
      };
      FIXED_VECTOR<3> hx[8];
      FIXED_VECTOR<NUM_COMPONENTS> hv[8];
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 8; k++) {
          hx[k] = x[hexi[j][k]];
          hv[k] = v[hexi[j][k]];
        }
        write_hexvue_hex_compressed (IOB, hx, hv);
      }
    }
  
 private: // object data ////////////////////////////////////////////

 private: // class members //////////////////////////////////////////
  
#ifdef min
# undef min
#endif
#ifdef max
# undef max
#endif
  double min (double a, double b) { if (a > b) return b; else return a; }
  double max (double a, double b) { if (a < b) return b; else return a; }

  int write_hexvue_hex_ascii (ckit_iobuf_t IOB,
                              POINT p[8],
                              FIXED_VECTOR<NUM_COMPONENTS> values[8]
                              )
    {
      int j, i, n = NUM_COMPONENTS;
      char buf[128];
      
#undef RUN_IOP
#define RUN_IOP(action)                                         \
   {                                                            \
     if (!action) {                                             \
        fprintf (stderr, "Error: %s\n", ((char *)ckit_io_error_explanation ()));     \
     }                                                          \
   }
      RUN_IOP (ckit_io_put_line (IOB, "HV"));
      FOR (j, 8) {
        sprintf (buf, "%g %g %g", p[j](0), p[j](1), p[j](2));
        RUN_IOP (ckit_io_put_line (IOB, buf));
        FOR (i, n) {
          sprintf (buf, "%g", values[j](i));
          RUN_IOP (ckit_io_put_line (IOB, buf));
        }      
      }
      
      return 1;
    }
  
  int write_hexvue_hex_compressed (ckit_iobuf_t IOB,
                                  POINT p[8],
                                  FIXED_VECTOR<NUM_COMPONENTS> values[8]
                                  )
    {
      int n = NUM_COMPONENTS;
      float flbox[6];
      unsigned char geombytes[24];
      float minv[NUM_COMPONENTS];
      float maxv[NUM_COMPONENTS];
      unsigned char byte[NUM_COMPONENTS];
      G3D_box_t box;
      
      for (int i = 0; i < n; i++) {
        minv[i] = FLT_MAX;
        maxv[i] = -FLT_MAX;
        for (int m = 0; m < 8; m++) {
          int j = m;
          minv[i] = min (minv[i], values[j][i]);
          maxv[i] = max (maxv[i], values[j][i]);
        }
      }
      
#define BX(WHAT, which)                                                         \
  {                                                                             \
    G3D_vect_t ap; ap.x = p[which](0); ap.y = p[which](1); ap.z = p[which](2);  \
    WHAT (box, ap);                                                             \
  }
      BX (G3D_INIT_BOX, 0);
      BX (G3D_UPDT_BOX, 1);
      BX (G3D_UPDT_BOX, 2);
      BX (G3D_UPDT_BOX, 3);
      BX (G3D_UPDT_BOX, 4);
      BX (G3D_UPDT_BOX, 5);
      BX (G3D_UPDT_BOX, 6);
      BX (G3D_UPDT_BOX, 7);
      
      flbox[0] = box.lower.x;
      flbox[1] = box.lower.y;
      flbox[2] = box.lower.z;
      flbox[3] = box.upper.x;
      flbox[4] = box.upper.y;
      flbox[5] = box.upper.z;
      
#undef RUN_IOP
#define RUN_IOP(action)                                         \
   {                                                            \
     if (!action) {                                             \
        fprintf (stderr, "Error: %s\n", ((char *)ckit_io_error_explanation ()));     \
     }                                                          \
   }
      RUN_IOP (ckit_io_put_line (IOB, "CHV"));
      RUN_IOP (ckit_io_put_data (IOB, minv,  n * sizeof (float)));
      RUN_IOP (ckit_io_put_data (IOB, maxv,  n * sizeof (float)));
      RUN_IOP (ckit_io_put_data (IOB, flbox, 6 * sizeof (float)));
      int k = 0;
      for (int m = 0; m < 8; m++) {
        int j = m;
        ckit_io_convert_double_to_byte (p[j](0), box.lower.x, box.upper.x,
                                        &geombytes[k]); k++;
        ckit_io_convert_double_to_byte (p[j](1), box.lower.y, box.upper.y,
                                        &geombytes[k]); k++;
        ckit_io_convert_double_to_byte (p[j](2), box.lower.z, box.upper.z,
                                        &geombytes[k]); k++;
      }
      RUN_IOP (ckit_io_put_data (IOB, geombytes, 24));
      for (int m = 0; m < 8; m++) {
        int j = m;
        for (int i = 0; i < n; i++) {
          ckit_io_convert_double_to_byte (values[j][i], minv[i], maxv[i],
                                          &byte[i]);
        }      
        RUN_IOP (ckit_io_put_data (IOB, byte, n * sizeof (unsigned char)));
      }
      
      return 1;
    }
  
};

template <int NUM_COMPONENTS>
static ECELL_HEXVUE<NUM_COMPONENTS> *
make (GCELL *gcell)
{
  return (new ECELL_HEXVUE_T10<NUM_COMPONENTS> (gcell));
}

template <int NUM_COMPONENTS>
bool
ECELL_HEXVUE_T10<NUM_COMPONENTS>::register_make_func ()
{
  return ECELL_HEXVUE<NUM_COMPONENTS>::register_make_func (make,
                                                           string (GCELL_SOLID_T10::TYPE_NAME),
                                                           string ("default"));
}

#endif
