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
#ifndef PROTO_XHEXVUE_H
# define PROTO_XHEXVUE_H

#include <list>
#include "algo.h"
#include "gmesh.h"
#include "ecell_xhexvue.h"
#include "proto_base.h"
#include "quant_evaluator.h"
extern "C" {
#include "ioP.h"
}

/**
This protocol writes data for all ecells 
that are associated with the protocol on input.
*/
class PROTO_XHEXVUE {

public: // object functions ////////////////////////////////////////

  
     PROTO_XHEXVUE (ALGO *algo){
      _algo = algo;      
     }
    ~PROTO_XHEXVUE () {}
    /**
    Writes data for all quantities whose names are given as 
    a list of strings for all ecells 
    that are associated with the protocol on input.
    */
    void write_xhexvue_data (ckit_iobuf_t IOB, PROTO_BASE *proto, list<string> quants) {
        SMART_HANDLE<PROTO_BASE::ECELL_ENUMERATOR> ee = proto->ecell_enumerator();
        ECELL *e;
        ee->reset();
        while ((e = ee->next())) {
          SMART_HANDLE<QUANT_EVALUATOR> qe = proto->quant_evaluator (e, quants);
          string path = "algorithms/xhexvue/" + _algo->name () + "/gcell_groups/" + e->gcell()->gcell_group()->name ();
          string implementation = _algo->mgr()->db()->DB_GET_STRING (path + "/implementation");
          ECELL_XHEXVUE *ecell_xhexvue = ECELL_XHEXVUE::make_ecell (e->gcell(), implementation);
          if (ecell_xhexvue){
            ecell_xhexvue->write (IOB, qe);
            delete(ecell_xhexvue);
          }
        }
    }   
 private:
    ALGO      *_algo; 
};



#endif
