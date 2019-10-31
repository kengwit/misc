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
#ifndef ALGO_HEXVUE_H
#  define ALGO_HEXVUE_H

#include <string>
#include "famexception.h"
#include "algo.h"
extern "C" {
#include "ioP.h"
}
#include "db.h"
#include "field.h"
#include "proto_hexvue.h"
#include "logger_stream_mgr.h"

template <int NUM_COMPONENTS>
class ALGO_HEXVUE : public ALGO {

 public: // object functions ////////////////////////////////////////

  /**
     Construct an algorithm.  Give it a name, and let it load
     its parameters from the database.
  */
  ALGO_HEXVUE (string name, MGR *mgr) : ALGO (name, mgr)
    {
      if (this->mgr()->db()) {
        string path = "algorithms/hexvue/" + name;
        _hexvue_file = this->mgr()->db()->DB_GET_STRING (path + "/hexvue_file");
      }
    }
  /**
     Run the algorithm.
  */
  void write_field (FIELD_VECTOR *display_on_geometry,
                    FIELD_VECTOR *eval_on_geometry,
                    FIELD<NUM_COMPONENTS> *field) {
    write_field_to_name (_hexvue_file, display_on_geometry, eval_on_geometry, field);
  }
  /**
     Run the algorithm.
  */
  void write_field_to_name (string filename,
                            FIELD_VECTOR *display_on_geometry,
                            FIELD_VECTOR *eval_on_geometry,
                            FIELD<NUM_COMPONENTS> *field);
  
 private: // object functions /////////////////////////////////////


 private: // object data //////////////////////////////////////////

  string _hexvue_file;
  
};

#include "logger_stream_op.h"
  /**
     Run the algorithm.
  */
template <int NUM_COMPONENTS>
  void ALGO_HEXVUE<NUM_COMPONENTS>::write_field_to_name (string filename,
                            FIELD_VECTOR *display_on_geometry,
                            FIELD_VECTOR *eval_on_geometry,
                            FIELD<NUM_COMPONENTS> *field) {
    char buf[1024];
    int n_output_values = 0;
    
    SMART_HANDLE<LOGGER_STREAM> lsb  = this->mgr()->logger("write_field_to_name", true);
    *lsb << "writing " << filename << endl;
    
    FILE *fp = ckit_io_open_stream ((char *) filename.c_str (), "w");
    CHECK (fp, EXCEPTION_NULL_PTR,;);
    
    ckit_iobuf_t IOB = ckit_io_new_buffer (fp);
    CHECK (IOB, EXCEPTION_NULL_PTR,;);
    
    ckit_io_put_line (IOB, "! HEXVUE file (C) 1995-2001 Petr Krysl.");
    ckit_io_put_line (IOB, "!   Generated by FAMULS (C) 2001 Petr Krysl.");
    sprintf (buf, "Field: %s", field->name ().c_str()); 
    ckit_io_put_line (IOB, buf);
    n_output_values += NUM_COMPONENTS; 
    
    sprintf (buf, "%d ! # of values per vertex", n_output_values);
    ckit_io_put_line (IOB, buf);
    
    for (int c = 0; c < NUM_COMPONENTS; c++) {
      sprintf (buf, "Field %s, component %d", (char *) field->name().c_str(), c);
      ckit_io_put_line (IOB, buf);
    }
    
    ckit_io_put_line (IOB, "! END OF HEADER");
    
    PROTO_HEXVUE<NUM_COMPONENTS> _proto_hexvue (this);
    _proto_hexvue.write_hexvue_data (IOB, display_on_geometry, eval_on_geometry, field);
    
    ckit_io_put_line (IOB, "END ! HEXVUE file.");
    
    ckit_io_delete_buffer (IOB);
    fclose (fp);
  }

#endif
