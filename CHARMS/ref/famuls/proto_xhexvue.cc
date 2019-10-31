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
#include "proto_xhexvue.h"

/*
// ashi: following function is already declared in h file. 
//if you have to open this comment, take care that pointer assignment is change. function defined in h file is corrected.
   
void 
PROTO_XHEXVUE::write_xhexvue_data (ckit_iobuf_t IOB, PROTO_BASE *proto, list<string> quants) 
{
    PROTO_BASE::ECELL_ENUMERATOR ee = proto->ecell_enumerator();
    ECELL *e;
    ee.reset();
    while (e = ee.next()) {
        QUANT_EVALUATOR qe = proto->quant_evaluator (e, quants);
        string path = "algorithms/xhexvue/" + _algo->name () + "/gcell_groups/" + e->gcell()->gcell_group->name ();
        string implementation = _algo->mgr()->db()->DB_GET_STRING (path + "/implementation");
        ECELL_XHEXVUE *ecell_xhexvue = ECELL_XHEXVUE::make_ecell (e->gcell(), implementation);
        ecell_xhexvue->write (IOB, qe);
    }
} 
  
*/
