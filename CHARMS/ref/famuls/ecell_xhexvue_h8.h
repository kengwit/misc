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
#ifndef ECELL_XHEXVUE_H8_H
# define ECELL_XHEXVUE_H8_H

#include "ecell_hexvue.h"
#include "gcell_solid_h8.h"
#include "ecell_xhexvue.h"

class ECELL_XHEXVUE_H8 : public ECELL_XHEXVUE {

 public: // declarations ////////////////////////////////////////////

  static bool register_make_func ();

 public: // class functions  ////////////////////////////////////////

  ECELL_XHEXVUE_H8 (GCELL *gcell)
    : ECELL_XHEXVUE::ECELL_XHEXVUE (gcell) {}
    ~ECELL_XHEXVUE_H8() {}
  /**
     Write the HEXVUE record for the ecell to the IO buffer.
   */
    // virtual void write (ckit_iobuf_t IOB, QUANT_EVALUATOR &qe);//ashi
    void write (ckit_iobuf_t IOB, SMART_HANDLE<QUANT_EVALUATOR>qe);
  
 private: // class members //////////////////////////////////////////
  
    int 
        write_hexvue_hex_ascii (ckit_iobuf_t IOB,
        POINT p[8],
        vector <double> vals[8]
        );
  
    int 
        write_hexvue_hex_compressed (ckit_iobuf_t IOB,
        POINT p[8],
        vector <double> vals[8]
        );

};


#endif
