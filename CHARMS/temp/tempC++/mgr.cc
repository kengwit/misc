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
#include "fen_maker.h"
#include "bfun_set.h"
#include "mgr.h"
#include "tokensP.h"
#include "expevalP.h"

int MGR::_instance_count =  0;

MGR::MGR (DB *db) 
{
  _db = db;
  CHECK (_db != 0, EXCEPTION_NULL_PTR,;);
  MGR::_instance_count++;
  CHECK (MGR::_instance_count == 1, EXCEPTION_ILLEGAL_USE,;);
  _mesh_mgr = new MESH_MGR (this);
  _mat_mgr  = new MAT_MGR (this);
  _load_mgr = new LOAD_MGR (this);

  // Define any expeval functions
  if (db->param_defined ("defines/funlist")) {
    list<string> tab_func_list = db->DB_GET_STRING_LIST ( "defines/funlist");
    tokens_parser parser = tokens_new_parser ();
    list<string>::iterator i = tab_func_list.begin();
    CHECK(i!= tab_func_list.end(),  EXCEPTION_NULL_PTR,;);
    string name;
    tokens_parse_line (parser, (char *) (*i).c_str ());
    if (tokens_total_of_tokens (parser) == 1) name =   tokens_token_as_string(parser,1);
    i++;
    CHECK(i!= tab_func_list.end(),  EXCEPTION_NULL_PTR,;);
    string file ;
    tokens_parse_line (parser, (char *) (*i).c_str ());
    if (tokens_total_of_tokens (parser) == 1) file =   tokens_token_as_string(parser,1);
    expeval_register_tabular_func ((char*) name.c_str(),(char*) file.c_str());
    tokens_delete_parser(parser);
  } 

}

MGR::~MGR () {
  if (_mesh_mgr) delete _mesh_mgr;
  if (_mat_mgr) delete _mat_mgr;
  if (_load_mgr) delete _load_mgr;
}

