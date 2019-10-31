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
#include "gcell_group.h"
#include "gsubmesh.h"
#include "gmesh.h"


GCELL_GROUP::GCELL_GROUP (string name, DB *db, GSUBMESH *gsubmesh)
{
  _name = name;
  _gsubmesh = gsubmesh;
  _gcell_type = "";
  // type
  {
    string path = "gcell_groups/" + _name + "/type";
    _gcell_type = db->DB_GET_STRING (path);
  }
  {
    string path = "gcell_groups/" + _name + "/conn_file";
    string gcell_file = db->DB_GET_STRING (path);
    CHECK (read_gcell_file (gcell_file), EXCEPTION_BAD_VALUE ,;);
  }
}

void 
GCELL_GROUP::add (GCELL *gcell)
{
  CHECK (_gcell_type == gcell->type_name (), EXCEPTION_BAD_VALUE,;);
  _gcells.push_back (gcell);
  gcell->set_gcell_group (this);
  //cerr << "adding "; gcell->debug_display ();
}

extern "C" {
#include "tokensP.h"  
}

bool
GCELL_GROUP::read_gcell_file (string file)
{
  tokens_parser parser = tokens_new_parser ();
  CHECK (parser, EXCEPTION_NULL_PTR,;);
  CHECK (tokens_open_file (parser, (char *) file.c_str ()), EXCEPTION_FILE_ERROR, (file));
  while (tokens_next_line (parser)) {
    int total_tok = tokens_total_of_tokens (parser);
    int nfens     = total_tok;
    if (tokens_have_keyword (parser, total_tok-1, "id")) nfens = total_tok - 2;
    vector <FEN *> fens;
    fens.resize (nfens);
    for (int i = 1; i <= nfens; i++) {
      int id = tokens_token_as_int (parser, i);
      fens[i-1] = _gsubmesh->gmesh()->find_fen (id);
      if (fens[i-1] == NULL) {
        cerr << "Node " << id << " not found" << endl;
        CHECK (fens[i-1], EXCEPTION_NULL_PTR,;);
      }
    }
    GCELL *gcell = GCELL::make_gcell (_gcell_type, "default", fens);
    if (!gcell) {
      cerr << "While reading \"" << file << "\": making " << _gcell_type << endl;
      CHECK (gcell, EXCEPTION_NULL_PTR,;);
    }
    if ( nfens == (total_tok-2 )) gcell->set_id (tokens_token_as_int (parser, total_tok));
    add (gcell);
  }
  tokens_close_curr_file (parser);
  tokens_delete_parser (parser);
  return true;
}

void
GCELL_GROUP::debug_display ()
{
  cerr << "GCELL_GROUP " << _name << endl;
  for (list <GCELL *>::const_iterator i = _gcells.begin (); i != _gcells.end (); i++) {
    GCELL *gcell = *i;
    gcell->debug_display ();
  }
}

GCELL_GROUP::~GCELL_GROUP () {
  for (list <GCELL *>::const_iterator i = _gcells.begin (); i != _gcells.end (); i++) {
    GCELL *gcell = *i;
    delete gcell;
  }
  _gcells.clear();
}
  
