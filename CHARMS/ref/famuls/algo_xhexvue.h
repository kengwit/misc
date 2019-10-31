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
#ifndef ALGO_XHEXVUE_H
#  define ALGO_XHEXVUE_H

#include <string>
#include "famexception.h"
#include "algo.h"
extern "C" {
#include "ioP.h"
#include "tokensP.h"
}
#include "db.h"
#include "field.h"
#include "proto_xhexvue.h"

class ALGO_XHEXVUE : public ALGO {

 public: // object functions ////////////////////////////////////////

  /**
     Construct an algorithm.  Give it a name, and let it load
     its parameters from the database.
  */
  

  ALGO_XHEXVUE (string name, MGR *mgr) : ALGO (name, mgr)
    {
      _name = name;
      tokens_parser parser = tokens_new_parser ();
      if (this->mgr()->db()) {
        string path = "algorithms/xhexvue/" + name;
        list <string> quants =  this->mgr()->db()->DB_GET_STRING_LIST (path + "/quants");
        
        for (list<string>::iterator ix = quants.begin(); ix != quants.end(); ix++) {
          tokens_parse_line (parser, (char *) (*ix).c_str ());
          string w  =   tokens_token_as_string(parser,1);
          _quants.push_back(w);
        }
      }
      tokens_delete_parser(parser);
    }
  /**
     Run the algorithm.
  */
  void write (PROTO_BASE *proto, double time);
  void write (PROTO_BASE *proto);

  void make_file_name(sizet cycle) {
    char t[50]; sprintf(t,"%ld",cycle);
    _xhexvue_file = _name + t + ".xvu";
  }

  void make_file_name(sizet cycle, string n) {
    char t[50]; sprintf(t,"%ld",cycle);
    _xhexvue_file = _name + t + n + ".xvu";
  }

  void make_file_name(string n) {
    _xhexvue_file = _name + n + ".xvu";
  }

 private: // object data //////////////////////////////////////////

  string        _xhexvue_file;
  list <string> _quants;
  string        _name; 
};

#endif
