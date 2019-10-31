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
#include "famexception.h"

using namespace std;

BASE_EXCEPTION::BASE_EXCEPTION () : _file(""), _line(0), _function(""), _condition(""), _exception_name("") {}


BASE_EXCEPTION::BASE_EXCEPTION (const char* f, const int l, const char *func,
                                const char* c, const char *e) :
  _file(f), _line(l), _function(func), _condition(c), _exception_name(e) {}


void BASE_EXCEPTION::set_fields (const char* f, const int l, const char *func, const char *c, const char *e)
{
  _file = f;
  _line = l;
  _function = func;
  _condition = c;
  _exception_name  = e;
}


void BASE_EXCEPTION::print_exception_data (ostream &out) const {
  out << "Error (" << _file << " l. " << _line << ", " << _function << "): " << endl;
  out << "  " << _condition << "; [" << _exception_name << "]" << endl;
}


void BASE_EXCEPTION::print_info (ostream &out) const {
  out << " " << endl;
}


