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
#ifndef ECELL_SILO_T10_H
# define ECELL_SILO_T10_H

#if defined (SILO) && SILO

#include "ecell_silo.h"
#include "gcell_solid_t10.h"

template <int NUM_COMPONENTS>
class ECELL_SILO_T10 : public ECELL_SILO<NUM_COMPONENTS> {

 public: // declarations ////////////////////////////////////////////

  static bool register_make_func ();

  static const SILO_WRITER::SILO_UCD_ZONE_SHAPE SHAPE = SILO_WRITER::ZONE_SHAPE_HEX;

 public: // class functions  ////////////////////////////////////////

  ECELL_SILO_T10 (GCELL *gcell)
    : ECELL_SILO<NUM_COMPONENTS>::ECELL_SILO (gcell) {}

  SILO_WRITER::SILO_UCD_ZONE_SHAPE shape () const { return SHAPE; }

  /**
     Write the SILO record for the ecell to the IO buffer.
     The cell will have the geometry given on input, and the field
     to display on the geometry is also specified.  To evaluate the
     field, the geometry eval_on_geometry is used.
   */
  virtual void stream_silo_data (SILO_WRITER *silo_writer,
                                 FIELD_VECTOR *display_on_geometry,
                                 FIELD_VECTOR *eval_on_geometry,
                                 FIELD<NUM_COMPONENTS> *field)
    {
      const int NFENS = gcell()->conn()->nfens();
      FIXED_VECTOR<3> x[15];
      FIXED_VECTOR<NUM_COMPONENTS> v[15];
      
      for (int pt = 0; pt < 15; pt++) {
        // evaluate the geometry
        POINT param_loc (g3dtet10_vertex_param_coords[pt][0],
                         g3dtet10_vertex_param_coords[pt][1],
                         g3dtet10_vertex_param_coords[pt][2]);
        EVALPT gevalpt (display_on_geometry, gcell(), param_loc);
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
        EVALPT fevalpt (field, gcell(), param_loc);
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

      silo_writer->start_zone (SILO_WRITER::ZONE_SHAPE_POLYHEDRON); {
        silo_writer->push_to_nodelist (4); // four faces
        for (sizet j = 0; j < 4; j++) { // for each face
          silo_writer->push_to_nodelist (6); // six vertices per face
          for (sizet i = 0; i < 6; i++) { // for each vertex
            sizet k = g3dtet10_vertex_on_face[j][i];
            SILO_WRITER::SILO_FEN_DATA fen_data;
            FEN *fen = gcell()->conn()->fen(k);
            fen_data.id = fen->id();
            fen_data.x  = x[k](0);
            fen_data.y  = x[k](1);
            fen_data.z  = x[k](2);
            for (int m = 0; m < NUM_COMPONENTS; m++) {
              fen_data.varvals[m] = v[k](m);
            }
            for (int m = NUM_COMPONENTS; m < SILO_WRITER::MAX_VARS; m++) {
              fen_data.varvals[m] = 0;
            }
            silo_writer->push_to_nodelist (&fen_data);
          } // for each vertex
        } // for each face
      } silo_writer->end_zone ();
    }
  
 private: // object data ////////////////////////////////////////////

};

template <int NUM_COMPONENTS>
static ECELL_SILO<NUM_COMPONENTS> *
make_silo_t10 (GCELL *gcell)
{
  return (new ECELL_SILO_T10<NUM_COMPONENTS> (gcell));
}

template <int NUM_COMPONENTS>
bool
ECELL_SILO_T10<NUM_COMPONENTS>::register_make_func ()
{
  return ECELL_SILO<NUM_COMPONENTS>::register_make_func (make_silo_t10,
                                                         string (GCELL_SOLID_T10::TYPE_NAME),
                                                         string ("default"));
}

#endif // SILO enabled?

#endif
