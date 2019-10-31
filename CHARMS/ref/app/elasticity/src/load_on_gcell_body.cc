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
#include "load_on_gcell_body.h"
extern "C" {
#include "tokensP.h"
}

const string LOAD_ON_GCELL_BODY::TYPE_NAME
         = LOAD_ON_GCELL::TYPE_NAME + " " +  "load_on_gcell_body";

static LOAD *
make (pair <DB *, string> p)
{
  return (new LOAD_ON_GCELL_BODY (p.first, p.second));
}

bool
LOAD_ON_GCELL_BODY::register_make_func ()
{
   tokens_parser _parser = tokens_new_parser ();
   const char *s = LOAD_ON_GCELL_BODY::TYPE_NAME.c_str();
   if (tokens_parse_line (_parser, (char *) s)) {
     for (int tok = 1; tok <= tokens_total_of_tokens (_parser); tok++) {
       //cerr << "Registering " << string(tokens_token_as_string (_parser, tok)) << endl;
       bool suc = LOAD::register_make_func (make, string(tokens_token_as_string (_parser, tok)));
       if (!suc) return false;
     }
   }
   tokens_delete_parser (_parser);
  return true;
}

LOAD_ON_GCELL_BODY::LOAD_ON_GCELL_BODY (DB *db, string name) : LOAD_ON_GCELL (db, name)
{
  string path = "loads/" + this->name();
  string expr[3];
  list<string> all = db->DB_GET_STRING_LIST (path + "/force_density");
  for (int j = 0; j < 3; j++) { expr[j] = ""; }
  int i = 0;
  for (list<string>::iterator j = all.begin(); j != all.end(); j++) {
    expr[i++] = *j;
  }
  _force_density = FUNC<3> (expr);
}

