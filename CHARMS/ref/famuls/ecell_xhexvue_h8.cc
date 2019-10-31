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
#include "ecell_xhexvue_h8.h"
#include "evalpt.h"
#include "fixed_vector.h"

static ECELL_XHEXVUE *
make (GCELL *gcell)
{
  return (new ECELL_XHEXVUE_H8 (gcell));
}

bool
ECELL_XHEXVUE_H8::register_make_func ()
{
    return ECELL_XHEXVUE::register_make_func (make,
        string (GCELL_SOLID_H8::TYPE_NAME),
        string ("default"));
}

int 
ECELL_XHEXVUE_H8::write_hexvue_hex_ascii (ckit_iobuf_t IOB,
                                  POINT p[8],
                                  vector <double> vals[8]
                                  )
  {
      const int n = vals[0].size();
      int j, i;
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
          sprintf (buf, "%g", vals[j][i]);
          RUN_IOP (ckit_io_put_line (IOB, buf));
        }      
      }
      
      return 1;
    }

int 
ECELL_XHEXVUE_H8::write_hexvue_hex_compressed (ckit_iobuf_t IOB,
                                  POINT p[8],
                                  vector <double> vals[8]
                                  )
    {
      const int n = vals[0].size();
      float flbox[6];
      unsigned char geombytes[24]; 
      float  *minv;
      float  *maxv;
      unsigned char  *byte;
      minv = new float[n];
      maxv = new float[n];
      byte = new unsigned char[n];

      G3D_box_t box;
      
      for (int i = 0; i < n; i++) {
        minv[i] = FLT_MAX;
        maxv[i] = -FLT_MAX;
        for (int m = 0; m < 8; m++) {
          int j = m;
          minv[i] = std::min ((double) minv[i], vals[j][i]);
          maxv[i] = std::max ((double) maxv[i], vals[j][i]);
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
#define RUN_IOP(action)                                                                 \
   {                                                                                    \
     if (!action) {                                                                     \
        fprintf (stderr, "Error: %s\n", ((char *)ckit_io_error_explanation ()));        \
     }                                                                                  \
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
          // ckit_io_convert_double_to_byte (values[j][i], minv[i], maxv[i], &byte[i]);//ashi : I might be wrong here.
          ckit_io_convert_double_to_byte (vals[j][i], minv[i], maxv[i], &byte[i]);
          
        }      
        RUN_IOP (ckit_io_put_data (IOB, byte, n * sizeof (unsigned char)));
      }
      delete [] minv; delete [] maxv;delete [] byte;
      return 1;
    }


void 
ECELL_XHEXVUE_H8::write (ckit_iobuf_t IOB, SMART_HANDLE<QUANT_EVALUATOR>qe) 
{
    const int NFENS = 8;
    FIXED_VECTOR<3> x[NFENS];
    vector <double> vals[NFENS];
    for (int pt = 0; pt < NFENS; pt++) {
        // evaluate the geometry
        POINT param_loc (g3dhex_vertex_param_coord[pt][0],
                         g3dhex_vertex_param_coord[pt][1],
                         g3dhex_vertex_param_coord[pt][2]);
        x[pt] = qe->loc(param_loc);
        vals[pt] = qe->quants(param_loc);
    } // for all pt
    write_hexvue_hex_ascii (IOB, x, vals);
}
  
