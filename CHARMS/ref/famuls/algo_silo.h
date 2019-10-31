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
#ifndef ALGO_SILO_H
#  define ALGO_SILO_H

#if defined (SILO) && SILO

#include <string>
#include "famexception.h"
#include "algo.h"
#include "silo_writer.h"
#include "db.h"
#include "field.h"
#include "proto_silo.h"

template <int NUM_COMPONENTS>
class ALGO_SILO : public ALGO {

 public: // object functions ////////////////////////////////////////

  /**
     Construct an algorithm.  Give it a name, and let it load
     its parameters from the database.
  */
  ALGO_SILO (string name, MGR *mgr) : ALGO (name, mgr)
    {
    }
  /**
     Run the algorithm.
  */
  void write_field_to_name (string filename,
                            FIELD_VECTOR *display_on_geometry,
                            FIELD_VECTOR *eval_on_geometry,
                            FIELD<NUM_COMPONENTS> *field) {
    char buf[1024];
    int n_output_values = 0;
    
    this->mgr()->logger() << "writing a SILO file: " << filename << endl;
    
    SILO_WRITER silo_writer;
    vector <string> var_names; var_names.clear();
    if (NUM_COMPONENTS == 1) {
      var_names.push_back (field->name());
    } else {
      for (int m = 0; m < NUM_COMPONENTS; m++) {
        char buf[6];
        sprintf (buf, "%d", m);
        var_names.push_back (field->name() + "_component_" + string(buf));
      }
    }
    silo_writer.setup (NUM_COMPONENTS, var_names);
    
    // Generate data
    PROTO_SILO<NUM_COMPONENTS> _proto_silo (this);
    _proto_silo.stream_silo_data (&silo_writer, display_on_geometry, eval_on_geometry, field);

    // Write it out
    silo_writer.write_silo (filename);
  }
  
 private: // object functions /////////////////////////////////////


 private: // object data //////////////////////////////////////////

};

#endif // Was SILO enabled?

#endif
