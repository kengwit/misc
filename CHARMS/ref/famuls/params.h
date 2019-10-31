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
#ifndef PARAMS_H
#  define PARAMS_H

#undef min
#undef max
#undef min
#include <string>
#include <fstream>
#include <iomanip>
#include <strstream>
#include <cstdlib>
#include <algorithm>
#include <list>
#include <map>
#include <vector>

#include "famexception.h"

using namespace std;

enum OUTPUT_FORMAT {
  TEXT_FORMAT, LaTeX_FORMAT, HTML_FORMAT
};

using namespace std;

struct TOKENS {
  
  static char comment_char;

  public:

  class TOKENS_BASE {
    
  public:
    
    virtual ~TOKENS_BASE () {};
    
    virtual bool is_equal (const string &test_string) const = 0;
    
    virtual string description () const = 0;
    
    virtual TOKENS_BASE *clone () const = 0;
    
  };
  
  class TOKENS_INTEGER : public TOKENS_BASE {
    
  public:
    
    virtual bool is_equal (const string &test_string) const;
    virtual string description () const;
    virtual TOKENS_BASE *clone () const;
    
  };
  
  class TOKENS_DOUBLE : public TOKENS_BASE {
    
  public:
    
    virtual bool is_equal (const string &test_string) const;
    virtual string description () const;
    virtual TOKENS_BASE *clone () const;
    
  };
  
  class TOKENS_CHOICE : public TOKENS_BASE {
    
  public:
    
    TOKENS_CHOICE (const string &seq);
    virtual bool is_equal (const string &test_string) const;
    virtual string description () const;
    virtual TOKENS_BASE *clone () const;
    
  private:
    
    string _value;
    
  };
  
  
  // Comma-separated list of values.  Note: commas are not allowed 
  // within the values themselves.
  class TOKENS_MULTIPLE_CHOICE : public TOKENS_BASE {
    
  public:
    
    TOKENS_MULTIPLE_CHOICE (const string &seq);
    virtual bool is_equal (const string &test_string) const;
    virtual string description () const;
    virtual TOKENS_BASE *clone () const;
    
  private:
    
    string _value;
    
  };
  
  class TOKENS_BOOL : public TOKENS_CHOICE {
    
  public:
    
    TOKENS_BOOL ();
    virtual TOKENS_BASE *clone () const;
    
  };
  
  class TOKENS_ANY : public TOKENS_BASE {
    
  public:
    
    TOKENS_ANY ();
    virtual bool is_equal (const string &test_string) const;
    virtual string description () const;
    virtual TOKENS_BASE *clone () const;
    
  };
  
};


class PARAMS_HANDLER {
  
 public:

  DECLARE_EXCEPTION_1_ARG (EXCEPTION_PARAMETER_NOT_FOUND, string, << arg1);

 public:
  
  PARAMS_HANDLER ();
  virtual ~PARAMS_HANDLER ();
  virtual bool read_input (istream &input);
  // Read lines from a char array.  The separator is the newline character ('\n').
  virtual bool read_input_from_string (const char *s);
  virtual bool read_input (const string &filename);
  bool no_errors() const;
  void clear ();
  bool dcl_param (const string &param,
                  const string &default_value,
                  const TOKENS::TOKENS_BASE &pattern = TOKENS::TOKENS_ANY());
  void enter_scope (const string &scope);
  bool leave_scope ();
  bool param_defined (const string &param_name);
  bool param_defined_from_path (const string &param_name);
  const string &get_string (const string &param_name) const;
  long int get_integer (const string &param_name) const;
  double get_double (const string &param_name) const;
  bool get_bool (const string &param_name) const;
  ostream & output_params (ostream &out, const OUTPUT_FORMAT format);
  void output_params_scope (ostream &out,
                            const OUTPUT_FORMAT format,
                            const unsigned int indent_level);
  list <string> get_scope_names () const;
  list <string> get_param_names () const;
  const string & get_string_from_path (const string &param_name);
  long int get_integer_from_path (const string &param_name);
  double get_double_from_path (const string &param_name);
  bool get_bool_from_path (const string &param_name);
  list <string> get_string_list (const string &param_name) const;
  list <string> get_string_list_from_path (const string &param_name);
  list <long int> get_integer_list (const string &param_name);
  list <long int> get_integer_list_from_path (const string &param_name);
  list <double> get_double_list (const string &param_name);
  list <double> get_double_list_from_path (const string &param_name);


  void set_scope_keyword (string scope_keyword) { _scope_keyword = scope_keyword; }
  
 private:
  
  bool _no_errors;
  string _scope_keyword;
  
  struct SCOPE {
    
    SCOPE (string name);

    SCOPE () {}

    ~SCOPE ();
    
    typedef map<string, pair<string,TOKENS::TOKENS_BASE*> > param_map_t;
    
    param_map_t _param_map;
    
    string _name;
    map<string, SCOPE*> _scopes;
    
  };
  
  vector<string> _scope_nesting;
  
  SCOPE _defaults;
  SCOPE _set;
  bool scan_line (string line, const unsigned int iline);
  void make_scope (const string &scope);
  SCOPE* get_curr_default_scope ();
  const SCOPE* get_curr_default_scope () const;
  SCOPE* get_curr_set_scope ();
  const SCOPE* get_curr_set_scope () const;
  void enter_scopes (list <string> l);
  void leave_scopes (list <string> l);

};


inline
bool
PARAMS_HANDLER::no_errors() const
{
  return _no_errors;
}


#endif

