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
#include "db.h"

DB::DB () {
  _params_handler = new PARAMS_HANDLER ();
  _params_handler->clear ();
  _params_handler->set_scope_keyword ("obj");
}

string DB::get_string_from_path (char *file, int line, const string param_path) {
  try {
    return _params_handler->get_string_from_path (param_path);
  } catch (EXCEPTION_BAD_VALUE e) {
    cerr << "Bad value: " << file << ", l." << line << " [" << param_path << "]" << endl;
    return 0;
  } catch (PARAMS_HANDLER::EXCEPTION_PARAMETER_NOT_FOUND) {
    cerr << "Parameter not found: " << file << ", l." << line << " [" << param_path << "]" << endl;
    return 0;
  }
}

double 
DB::get_double_from_path (char *file, int line, const string param_path) {
  try{
    return _params_handler->get_double_from_path (param_path);
  } catch (EXCEPTION_BAD_VALUE e) {
    cerr << "Bad value: " << file << ", l." << line << " [" << param_path << "]" << endl;
    return 0;
  } catch (PARAMS_HANDLER::EXCEPTION_PARAMETER_NOT_FOUND) {
    cerr << "Parameter not found: " << file << ", l." << line << " [" << param_path << "]" << endl;
    return 0;
  }
}

long int 
DB::get_integer_from_path (char *file, int line, const string param_path) {
  try{
    return _params_handler->get_integer_from_path (param_path);
  } catch (EXCEPTION_BAD_VALUE e) {
    cerr << "Bad value: " << file << ", l." << line << " [" << param_path << "]" << endl;
    return 0;
  } catch (PARAMS_HANDLER::EXCEPTION_PARAMETER_NOT_FOUND) {
    cerr << "Parameter not found: " << file << ", l." << line << " [" << param_path << "]" << endl;
    return 0;
  }
}

bool 
DB::get_bool_from_path (char *file, int line, const string param_path) {
  try{
    return _params_handler->get_bool_from_path (param_path);
  } catch (EXCEPTION_BAD_VALUE e) {
    cerr << "Bad value: " << file << ", l." << line << " [" << param_path << "]" << endl;
    return 0;
  } catch (PARAMS_HANDLER::EXCEPTION_PARAMETER_NOT_FOUND) {
    cerr << "Parameter not found: " << file << ", l." << line << " [" << param_path << "]" << endl;
    return 0;
  }
}

bool
DB::load_file (string file)
{
  ifstream ifs (file.c_str());
  return _params_handler->read_input (ifs);
}

bool
DB::load_string (string s)
{
  return _params_handler->read_input_from_string (s.c_str());
}

list <string>
DB::get_string_list_from_path (char *file, int line, const string param_path)
{
  try{
    return _params_handler->get_string_list_from_path (param_path);
  } catch (EXCEPTION_BAD_VALUE e) {
    cerr << "Bad value: " << file << ", l." << line << " [" << param_path << "]" << endl;
    list <string> l; l.clear ();
    return l;
  } catch (PARAMS_HANDLER::EXCEPTION_PARAMETER_NOT_FOUND) {
    cerr << "Parameter not found: " << file << ", l." << line << " [" << param_path << "]" << endl;
    list <string> l;
    return l;
  }
}

list <long int>
DB::get_integer_list_from_path (char *file, int line, const string param_path)
{
  try{
    return _params_handler->get_integer_list_from_path (param_path);
  } catch (EXCEPTION_BAD_VALUE e) {
    cerr << "Bad value: " << file << ", l." << line << " [" << param_path << "]" << endl;
    list <long int> l; l.clear ();
    return l;
  } catch (PARAMS_HANDLER::EXCEPTION_PARAMETER_NOT_FOUND) {
    cerr << "Parameter not found: " << file << ", l." << line << " [" << param_path << "]" << endl;
    list <long int> l;
    return l;
  }
}

list <double>
DB::get_double_list_from_path (char *file, int line, const string param_path)
{
  try{
    return _params_handler->get_double_list_from_path (param_path);
  } catch (EXCEPTION_BAD_VALUE e) {
    cerr << "Bad value: " << file << ", l." << line << " [" << param_path << "]" << endl;
    list <double> l; l.clear ();
    return l;
  } catch (PARAMS_HANDLER::EXCEPTION_PARAMETER_NOT_FOUND) {
    cerr << "Parameter not found: " << file << ", l." << line << " [" << param_path << "]" << endl;
    list <double> l;
    return l;
  }
}

bool
DB::param_defined (const string param_path)
{
  return _params_handler->param_defined_from_path (param_path);
}

DB::~DB () {
  if (_params_handler) delete _params_handler;
}
