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
#ifndef FENSRC_H
# define FENSRC_H

#include "fen.h"
#include "conn.h"
#include "conn_ref.h"
#include <list>
#include <map>
#include "fen_maker.h"
#include "limits.h"

class REF_FEN_SRC {

 public: // object methods //////////////////////////////////////////

  REF_FEN_SRC ();
  ~REF_FEN_SRC () { // RAW need to cleanup the map & the conn_ref's
    for (conn_ref_list_map_map_t::iterator i = _crlmm.begin();
         i != _crlmm.end(); i++) {
      conn_ref_list_map_t *crlm = i->second;
      for (conn_ref_list_map_t::iterator j = crlm->begin(); j != crlm->end(); j++) {
        conn_ref_list_t *l = j->second;
        for (conn_ref_list_t::iterator k = l->begin(); k != l->end(); k++) {
          CONN_REF *conn_ref = *k;
          delete conn_ref;
        }
        delete l;
      }
      delete crlm;
    }
  }
  /**
     Get a refinement node.  All connectivities defined so far
     are searched for a match.  If the match is found, the refinement
     node index indx is mapped into an index corresponding to the
     re-oriented connectivity, and the refinement node is returned;
     otherwise the connectivity
     is defined to be recognized by any subsequent searches, a new refinement
     node is generated and returned.
  */
  FEN * get_ref_fen (CONN_BASE *conn, sizet indx, FEN_MAKER *fen_maker);
  
 private: // object data ////////////////////////////////////////////

  // This collects all conn_ref objects tagged with the same node number.
  typedef list <CONN_REF *> conn_ref_list_t;
  // This maps the lowest node number in the connectivity to the conn_ref_list.
  typedef map <sizet, conn_ref_list_t *> conn_ref_list_map_t;
  // This maps the conn_code to the 
  typedef map <CONN_CODE, conn_ref_list_map_t *> conn_ref_list_map_map_t;

  conn_ref_list_map_map_t _crlmm;

 private: // object methods /////////////////////////////////////////

  CONN_REF *find_or_make_conn_ref (CONN_BASE *conn);

};

#endif
