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
#include "params.h"


char TOKENS::comment_char = '!';


bool
TOKENS::TOKENS_INTEGER::is_equal (const string &test_string) const
{
  istrstream str(test_string.c_str());
  int i;
  if (str >> i) return true;
  return false;
}


string
TOKENS::TOKENS_INTEGER::description () const {return "integer";}


TOKENS::TOKENS_BASE *
TOKENS::TOKENS_INTEGER::clone () const {return new TOKENS::TOKENS_INTEGER();}


bool
TOKENS::TOKENS_DOUBLE::is_equal (const string &test_string) const
{
  istrstream str(test_string.c_str());
  double d;
  if (str >> d) return true;
  return false;
}


string
TOKENS::TOKENS_DOUBLE::description () const {return "double";}


TOKENS::TOKENS_BASE *
TOKENS::TOKENS_DOUBLE::clone () const {return new TOKENS::TOKENS_DOUBLE();}


TOKENS::TOKENS_CHOICE::TOKENS_CHOICE (const string &v)
{
  _value = v;
  while (_value.find(" |") != string::npos) _value.replace (_value.find(" |"), 2, "|");
  while (_value.find("| ") != string::npos) _value.replace (_value.find("| "), 2, "|");
}


bool
TOKENS::TOKENS_CHOICE::is_equal (const string &test_string) const
{
  vector<string> choices;
  string tmp(_value);
  while (tmp.find('|') != string::npos) {
    if (test_string == string(tmp, 0, tmp.find('|')))
      return true;
    tmp.erase (0, tmp.find('|')+1);
  }
  if (test_string == tmp)
    return true;
  return false;
}


string
TOKENS::TOKENS_CHOICE::description () const {return _value;}


TOKENS::TOKENS_BASE *
TOKENS::TOKENS_CHOICE::clone () const {return new TOKENS::TOKENS_CHOICE(_value);}


TOKENS::TOKENS_MULTIPLE_CHOICE::TOKENS_MULTIPLE_CHOICE (const string &val)
{
  CHECK (val.find (",") == string::npos, EXCEPTION_BAD_VALUE,;);
  
  _value = val;
  while (_value.find(" |") != string::npos)
    _value.replace (_value.find(" |"), 2, "|");
  while (_value.find("| ") != string::npos)
    _value.replace (_value.find("| "), 2, "|");
}


bool
TOKENS::TOKENS_MULTIPLE_CHOICE::is_equal (const string &test_string_list) const {
  string tmp = test_string_list;
  list<string> split_list;
  
  // first split the input list
  while (tmp.length() != 0) {
    string name;
    name = tmp;
    
    if (name.find(",") != string::npos) {
      name.erase (name.find(","), string::npos);
      tmp.erase (0, test_string_list.find(",")+1);
    } else
      tmp = "";
    
    while ((name.length() != 0) &&
           (name[0] == ' '))
      name.erase (0,1);
    while (name[name.length()-1] == ' ')
      name.erase (name.length()-1, 1);
    
    split_list.push_back (name);
  }
  
  for (list<string>::const_iterator test_string = split_list.begin();
       test_string != split_list.end(); ++test_string) {
    bool string_found = false;
    tmp = _value;
    while (tmp.find('|') != string::npos) {
      if (*test_string == string(tmp, 0, tmp.find('|'))) {
        string_found = true;
        break;
      }
      tmp.erase (0, tmp.find('|')+1);
    }
    if (!string_found)
      if (*test_string == tmp)
        string_found = true;
    if (!string_found)
      return false;
  }
  
  return true;
}


string
TOKENS::TOKENS_MULTIPLE_CHOICE::description () const {return _value;}


TOKENS::TOKENS_BASE *
TOKENS::TOKENS_MULTIPLE_CHOICE::clone () const {return new TOKENS::TOKENS_MULTIPLE_CHOICE(_value);}


TOKENS::TOKENS_BOOL::TOKENS_BOOL () : TOKENS_CHOICE ("true|false") {}


TOKENS::TOKENS_BASE *
TOKENS::TOKENS_BOOL::clone () const {return new TOKENS::TOKENS_BOOL();}


TOKENS::TOKENS_ANY::TOKENS_ANY () {}


bool
TOKENS::TOKENS_ANY::is_equal (const string &) const {return true;}


string
TOKENS::TOKENS_ANY::description () const {return "any";}


TOKENS::TOKENS_BASE *
TOKENS::TOKENS_ANY::clone () const {return new TOKENS::TOKENS_ANY();}


PARAMS_HANDLER::PARAMS_HANDLER () : _no_errors(true),  _scope_keyword("scope") {}


PARAMS_HANDLER::~PARAMS_HANDLER () {}


bool
PARAMS_HANDLER::read_input (istream &input)
{
  CHECK (input != 0, EXCEPTION_NULL_PTR,;);

  string line;
  int linenum=0;
  while (input) {
    ++linenum;
    getline (input, line);
    if (!scan_line (line, linenum)) 
      _no_errors = false;
  }
  
  return _no_errors;
}


bool
PARAMS_HANDLER::read_input (const string &filename)
{
  ifstream input (filename.c_str());

  if (!input) {
    cerr << "PARAMS_HANDLER::read_input: could not open file <"
         << filename << "> for reading." << endl
         << "Trying to make file \""
         << filename << "\" with default values for you." << endl;
    
    ofstream output (filename.c_str());
    if (output) {
      output_params (output, TEXT_FORMAT);
    }
    
    return false;
  }
  
  return read_input (input);
}


bool
PARAMS_HANDLER::read_input_from_string (const char *s)
{
  if ((s == 0) || ((*s) == 0)) return true;
  
  string line;
  string input (s);
  int    linenum=0;
  
  if (input[input.length()-1] != '\n')
    input += '\n';
  
  while (input.size() != 0) {
    line.assign (input, 0, input.find('\n'));
    input.erase (0, input.find('\n')+1);
    ++linenum;
    if (!scan_line (line, linenum)) 
      _no_errors = false;
  }
  
  return _no_errors;
}


void
PARAMS_HANDLER::clear ()
{
  _no_errors = true;
  
  _scope_nesting.clear ();
  _defaults._param_map.clear ();
  _set._param_map.clear ();
  
  map<string, SCOPE*>::iterator p;
  
  for (p = _defaults._scopes.begin(); p != _defaults._scopes.end(); ++p)
    delete p->second;
  for (p = _set._scopes.begin(); p != _set._scopes.end(); ++p) {
    if (p->second) {
      delete p->second;
      CHECK (false, EXCEPTION_BAD_VALUE,;);
    }
  }
  _defaults._scopes.clear ();
  _set._scopes.clear ();
}


bool
PARAMS_HANDLER::dcl_param (const string &param_name,
                                const string &default_value,
                                const TOKENS::TOKENS_BASE &pattern)
{
  SCOPE *p = get_curr_default_scope ();
  
  CHECK (p->_param_map.find (param_name) == p->_param_map.end(), EXCEPTION_BAD_VALUE,;);

  CHECK (pattern.is_equal (default_value), EXCEPTION_BAD_VALUE,;);
  
  if (p->_param_map.find (param_name) != p->_param_map.end())
    return false;
  
  p->_param_map[param_name] = make_pair(default_value, pattern.clone());
  
  if (!pattern.is_equal (default_value))
    return false;
  
  return true;
}


void
PARAMS_HANDLER::make_scope (const string &scope_name) 
{
  SCOPE *pd = get_curr_default_scope ();
  
  if (pd->_scopes.find (scope_name) == pd->_scopes.end()) {
    pd->_scopes[scope_name] = new SCOPE(scope_name);
    
    SCOPE *pc = get_curr_set_scope ();
    pc->_scopes[scope_name] = new SCOPE(scope_name);
  }
}

void
PARAMS_HANDLER::enter_scope (const string &scope) 
{
  make_scope (scope);
  _scope_nesting.push_back (scope);
}


bool
PARAMS_HANDLER::leave_scope ()
{
  CHECK (_scope_nesting.size() != 0, EXCEPTION_BAD_SYNTAX,;);
  if (_scope_nesting.size() == 0) 
    return false;
  _scope_nesting.pop_back ();
  return true;
}

typedef struct {
  list <string> path;
  string        param_name;
}  param_path_t;

static void
param_path (string param_path,  param_path_t *pp)
{
  pp->path.clear ();
  string tmp = string (param_path);
  sizet slashindx = 0;
  if ((slashindx = tmp.find ("/")) == 0) tmp.replace (slashindx, slashindx+1, "");
  while ((slashindx = tmp.find ("/")) < tmp.size ()) {
    string sn = string (tmp, 0, slashindx);
    pp->path.push_back (sn);
    tmp.replace (0, slashindx+1, "");
  }
  pp->param_name = tmp;
}

void
PARAMS_HANDLER::enter_scopes (list <string> l)
{
  list <string>::const_iterator i = l.begin ();
  while (i != l.end ()) {
    enter_scope ((*i));
    i++;
  }
}

void
PARAMS_HANDLER::leave_scopes (list <string> l)
{
  list <string>::const_iterator i = l.begin ();
  while (i != l.end ()) {
    leave_scope ();
    i++;
  }
}

const string &
PARAMS_HANDLER::get_string_from_path (const string &param_name) 
{
  param_path_t pp;
  param_path (param_name, &pp);
  enter_scopes (pp.path);
  try {
    const string &result = get_string (pp.param_name);
    leave_scopes (pp.path);
    return result;
  } catch (...) {
    leave_scopes (pp.path);
    throw;
  }
}

const string &
PARAMS_HANDLER::get_string (const string &param_name) const
{
  const PARAMS_HANDLER::SCOPE *pd = get_curr_default_scope ();
  const PARAMS_HANDLER::SCOPE *pc = get_curr_set_scope ();
  
  CHECK_THROW (pd->_param_map.find (param_name) != pd->_param_map.end(),
               EXCEPTION_PARAMETER_NOT_FOUND, ("get_string"));
  
  SCOPE::param_map_t::const_iterator ptr;
  ptr = pc->_param_map.find (param_name);
  if (ptr != pc->_param_map.end())
    return ptr->second.first;
  
  ptr = pd->_param_map.find (param_name);
  return ptr->second.first;
}

bool
PARAMS_HANDLER::param_defined (const string &param_name) 
{
  const PARAMS_HANDLER::SCOPE *pd = get_curr_default_scope ();
  
  return (pd->_param_map.find (param_name) != pd->_param_map.end());
}

bool
PARAMS_HANDLER::param_defined_from_path (const string &param_name) 
{
  param_path_t pp;
  param_path (param_name, &pp);
  enter_scopes (pp.path);
    bool result = param_defined (pp.param_name);
  leave_scopes (pp.path);
  return result;
}

static list <string>
value_list (string param_value)
{
  list <string> l;
  l.clear ();
  string tmp = string (param_value);
  sizet slashindx = 0;
  while ((slashindx = tmp.find (" ")) < tmp.size ()) {
    string sn = string (tmp, 0, slashindx);
    l.push_back (sn);
    tmp.replace (0, slashindx+1, "");
  }
  l.push_back (tmp);
  return l;
}

list <string>
PARAMS_HANDLER::get_string_list (const string &param_name) const
{
  try {
    string s = get_string (param_name);
    list <string> l = value_list (s);
    return l;
  } catch (EXCEPTION_PARAMETER_NOT_FOUND) {
    throw;
  }
}

list <string>
PARAMS_HANDLER::get_string_list_from_path (const string &param_name) 
{
  list <string> result;
  param_path_t pp;
  param_path (param_name, &pp);
  enter_scopes (pp.path); 
  try {
    result = get_string_list (pp.param_name);
    leave_scopes (pp.path);
    return result;
  }  catch (EXCEPTION_PARAMETER_NOT_FOUND) {
    leave_scopes (pp.path);
    throw;
  }
}

long int
PARAMS_HANDLER::get_integer_from_path (const string &param_name) 
{
  long int result;
  param_path_t pp;
  param_path (param_name, &pp);
  enter_scopes (pp.path);
    result = get_integer (pp.param_name);
  leave_scopes (pp.path);
  return result;
}


long int
PARAMS_HANDLER::get_integer (const string &param_name) const
{
  try {
    string s = get_string (param_name);
    char *endptr;
    long int i = strtol (s.c_str(), &endptr, 10);
    CHECK_THROW ((s.c_str()!='\0') || (*endptr == '\0'), EXCEPTION_BAD_VALUE,;);
    return i;
  } catch (EXCEPTION_PARAMETER_NOT_FOUND) {
    throw;
  } 
}

list <long int>
PARAMS_HANDLER::get_integer_list (const string &param_name) 
{
  try {
    string s = get_string (param_name);
    list <string> l = value_list (param_name);
    list <long int> result;
    list <string>::const_iterator ii = l.begin ();
    while (ii != l.end ()) {
      char *endptr;
      string s = *ii;
      long int i = strtol (s.c_str(), &endptr, 10);
      CHECK_THROW ((s.c_str()!='\0') || (*endptr == '\0'), EXCEPTION_BAD_VALUE,;);
      result.push_back (i);
      ii++;
    }
    return result;
  }  catch (EXCEPTION_PARAMETER_NOT_FOUND) {
    throw;
  }
}

list <long int>
PARAMS_HANDLER::get_integer_list_from_path (const string &param_path) 
{
  try {
    string s = get_string_from_path (param_path);
    list <string> l = value_list (s);
    list <long int> result;
    list <string>::const_iterator ii = l.begin ();
    while (ii != l.end ()) {
      char *endptr;
      string s = *ii;
      long int i = strtol (s.c_str(), &endptr, 10);
      CHECK_THROW ((s.c_str()!='\0') || (*endptr == '\0'), EXCEPTION_BAD_VALUE,;);
      result.push_back (i);
      ii++;
    }
    return result;
  } catch (EXCEPTION_PARAMETER_NOT_FOUND) {
    throw;
  }
}


double
PARAMS_HANDLER::get_double_from_path (const string &param_name) 
{
  double result;
  param_path_t pp;
  param_path (param_name, &pp);
  enter_scopes (pp.path);
    result = get_double (pp.param_name);
  leave_scopes (pp.path);
  return result;
}

list <double>
PARAMS_HANDLER::get_double_list (const string &param_name) 
{
  try {
    string s = get_string (param_name);
    list <string> l = value_list (s);
    list <double> result;
    list <string>::const_iterator ii = l.begin ();
    while (ii != l.end ()) {
      char *endptr;
      string s = *ii;
      double d = strtod (s.c_str(), &endptr);
      CHECK_THROW ((s.c_str()!='\0') || (*endptr == '\0'), EXCEPTION_BAD_VALUE,;);
      result.push_back (d);
      ii++;
    }
    return result;
  } catch (EXCEPTION_PARAMETER_NOT_FOUND) {
    throw;
  }
}

list <double>
PARAMS_HANDLER::get_double_list_from_path (const string &param_path) 
{
  try {
    string s = get_string_from_path (param_path);
    list <string> l = value_list (s);
    list <double> result;
    list <string>::const_iterator ii = l.begin ();
    while (ii != l.end ()) {
      char *endptr;
      string s = *ii;
      double d = strtod (s.c_str(), &endptr);
      CHECK_THROW ((s.c_str()!='\0') || (*endptr == '\0'), EXCEPTION_BAD_VALUE,;);
      result.push_back (d);
      ii++;
    }
    return result;
  } catch (EXCEPTION_PARAMETER_NOT_FOUND) {
    throw;
  }
}

double
PARAMS_HANDLER::get_double (const string &param_name) const
{
  try {
    string s = get_string (param_name);
    char *endptr;
    double d = strtod (s.c_str(), &endptr);
    CHECK_THROW ((s.c_str()!='\0') || (*endptr == '\0'), EXCEPTION_BAD_VALUE,;);
    
    return d;
  }  catch (EXCEPTION_PARAMETER_NOT_FOUND) {
    throw;
  }
}


bool
PARAMS_HANDLER::get_bool (const string &param_name) const
{
  try {
    string s = get_string (param_name);
    
    CHECK_THROW ((s=="true") || (s=="false"), EXCEPTION_BAD_VALUE,;);
    if (s=="true")
      return true;
    else
      return false;
  }  catch (EXCEPTION_PARAMETER_NOT_FOUND) {
    throw;
  }
}


bool
PARAMS_HANDLER::get_bool_from_path (const string &param_name) 
{
  try {
    bool result;
    param_path_t pp;
    param_path (param_name, &pp);
    enter_scopes (pp.path);
    result = get_bool (pp.param_name);
    leave_scopes (pp.path);
    return result;
  } catch (EXCEPTION_PARAMETER_NOT_FOUND) {
    throw;
  }
}


ostream &
PARAMS_HANDLER::output_params (ostream &out, OUTPUT_FORMAT format)
{
  CHECK ((format == TEXT_FORMAT) || (format == LaTeX_FORMAT), EXCEPTION_NOT_IMPLEMENTED,;);
  
  CHECK (out != 0, EXCEPTION_NULL_PTR,;);
  
  switch (format) {
  case TEXT_FORMAT:
    out << "!Parameters:" << endl
        << "!-----------" << endl;
    break;
  case LaTeX_FORMAT:
    out << "\\subsubsection*{Parameters}";
    out << endl << endl;
    out << "\\begin{itemize}"
        << endl;
    break;
  default:
    CHECK (false, EXCEPTION_NOT_IMPLEMENTED,;);
  }
  
  output_params_scope (out, format, 0);
  
  switch (format) {
    case TEXT_FORMAT:
      break;
  case LaTeX_FORMAT:
    out << "\\end{itemize}" << endl;
    break;
  default:
    CHECK (false, EXCEPTION_NOT_IMPLEMENTED,;);
  }
  
  return out;
}


void
PARAMS_HANDLER::output_params_scope (ostream           &out,
                                          const OUTPUT_FORMAT  format,
                                          const unsigned int indent_level) 
{
  CHECK ((format == TEXT_FORMAT) || (format == LaTeX_FORMAT), EXCEPTION_NOT_IMPLEMENTED,;);
  
  CHECK (out != 0, EXCEPTION_NULL_PTR,;);
  
  SCOPE *pd = get_curr_default_scope ();
  SCOPE *pc = get_curr_set_scope ();
  
  // traverse param list
  SCOPE::param_map_t::const_iterator ptr;
  
  // first find out the longest param name
  unsigned int longest_param_name = 0;
  for (ptr = pd->_param_map.begin(); ptr != pd->_param_map.end(); ++ptr)
    if (ptr->first.length() > longest_param_name)
      longest_param_name = ptr->first.length();
  
  for (ptr = pd->_param_map.begin(); ptr != pd->_param_map.end(); ++ptr) {
    if ((pc->_param_map.find(ptr->first) != pc->_param_map.end()) &&
        (pc->_param_map[ptr->first].first != pd->_param_map[ptr->first].first))
      switch (format) {
      case TEXT_FORMAT:
        out << setw(indent_level*2) << ""
            << "param "
            << ptr->first
            << setw(longest_param_name-ptr->first.length()+3) << " = "
            << pc->_param_map[ptr->first].first
            << "  ! "
            << pd->_param_map[ptr->first].first
            << endl;
        break;
      case LaTeX_FORMAT:
        out << "\\item {\\bf " << ptr->first << ":} "
            << pc->_param_map[ptr->first].first
            << " ({\\it default:} "
            << pd->_param_map[ptr->first].first
            << ")"
            << endl;
        break;
      default:
        CHECK (false, EXCEPTION_NOT_IMPLEMENTED,;);
      }
    // not a changed param
    else
      switch (format) {
      case TEXT_FORMAT:
        out << setw(indent_level*2) << ""
            << "param "
            << ptr->first
            << setw(longest_param_name-ptr->first.length()+3) << "= "
            << ptr->second.first << endl;
        break;
      case LaTeX_FORMAT:
        out << "\\item {\\bf " << ptr->first << ":} "
            << ptr->second.first
            << endl;
        break;
      default:
        CHECK (false, EXCEPTION_NOT_IMPLEMENTED,;);
      }
  }
  
  
  map<string, SCOPE*>::const_iterator ptrss;
  for (ptrss = pd->_scopes.begin(); ptrss != pd->_scopes.end(); ++ptrss)
    {
      switch (format) 
	{
        case TEXT_FORMAT:
          out << setw(indent_level*2) << ""
              << _scope_keyword << " " << ptrss->first << endl;
          break;
        case LaTeX_FORMAT:
          out << endl
              << "\\item {\\bf "
              << _scope_keyword << " " << ptrss->first
              << "}" << endl
              << "\\begin{itemize}"
              << endl;
          break;
        default:
          CHECK (false, EXCEPTION_NOT_IMPLEMENTED,;);
	}
      enter_scope (ptrss->first);
      output_params_scope (out, format, indent_level+1);
      leave_scope ();
      switch (format) {
      case TEXT_FORMAT:
        out << setw(indent_level*2) << ""
            << "end" << endl
            << endl;
        if (indent_level == 0)
          out << endl;
        break;
      case LaTeX_FORMAT:
        out << "\\end{itemize}"
            << endl;
        break;
      default:
        CHECK (false, EXCEPTION_NOT_IMPLEMENTED,;);
        break;
      }
    }
}


bool
PARAMS_HANDLER::scan_line (string line, const unsigned int linenum)
{
  if (line.find(TOKENS::comment_char) != string::npos)
    line.erase (line.find(TOKENS::comment_char), string::npos);
  while (line.find('\t') != string::npos)
    line.replace (line.find('\t'), 1, " ");
  while (line.find("  ") != string::npos)
    line.erase (line.find("  "), 1);
  if ((line.length() != 0) && (line[0] == ' '))  line.erase (0, 1);
  // empty line? quit
  if (line.length() == 0) return true;
  
  if (line[line.length()-1] == ' ')
    line.erase (line.size()-1, 1);

  if ((line.find (_scope_keyword + " ") == 0)) {
    line.erase (0, _scope_keyword.length()+1);
    
    string scopename = line;
    SCOPE *pc = get_curr_set_scope ();
    if (pc->_scopes.find(scopename) == pc->_scopes.end()) {
      make_scope (scopename);
    }
    
    _scope_nesting.push_back (scopename);
    return true;
  }
  
  // leave scope
  if ((line.find ("END") == 0) ||
      (line.find ("End") == 0) ||
      (line.find ("end") == 0))
    if (_scope_nesting.size() == 0) {
      cerr << "Line " << linenum
           << ": No scope to leave." << endl;
      return false;
    } else
      return leave_scope ();
  
  // param statement
  if ((line.find ("PARAM ") == 0) ||
      (line.find ("Param ") == 0) ||
      (line.find ("param ") == 0)) {
    // Erase the PARAM statement and eliminate
    // spaces around the '='
    line.erase (0, 6);
    if (line.find(" =") != string::npos)
      line.replace (line.find(" ="), 2, "=");
    if (line.find("= ") != string::npos)
      line.replace (line.find("= "), 2, "=");
    
    // extract param name and value
    string param_name  (line, 0, line.find('='));
    string param_value (line, line.find('=')+1, string::npos);
    
    SCOPE *pd = get_curr_default_scope ();
    
    if (pd->_param_map.find(param_name) == pd->_param_map.end()) {
      dcl_param (param_name, "(no default value)", TOKENS::TOKENS_ANY());
    }
    
    if (param_value.find ('{') == string::npos)
      if (!pd->_param_map[param_name].second->is_equal(param_value)) {
          cerr << "Line " << linenum << ":" << endl
               << "    Value" << endl
               << "        " << param_value << endl
               << "    for the param" << endl
               << "        " << param_name << endl
               << "    is not valid: " << endl
               << "        " << pd->_param_map[param_name].second->description() << endl;
          return false;
        }
    
      SCOPE *pc = get_curr_set_scope ();
      pc->_param_map[param_name] = make_pair(param_value,
					  static_cast<TOKENS::TOKENS_BASE*>(0));
      
      return true;
    }
  
  // this line matched nothing known
  cerr << "Line " << linenum << ": Not recognized:" << endl
       << "    " << line << endl;
  return false;
}


PARAMS_HANDLER::SCOPE *
PARAMS_HANDLER::get_curr_default_scope () 
{
  SCOPE *scope = &_defaults;
  vector<string>::const_iterator scopename = _scope_nesting.begin();
  
  while (scopename != _scope_nesting.end()) {
    scope = scope->_scopes[*scopename];
    ++scopename;
  }
  
  return scope;
}


const PARAMS_HANDLER::SCOPE *
PARAMS_HANDLER::get_curr_default_scope () const 
{
  SCOPE *scope = const_cast<SCOPE*>(&_defaults);
  
  vector<string>::const_iterator scopename = _scope_nesting.begin();
  
  while (scopename != _scope_nesting.end()) {
    scope = scope->_scopes[*scopename];
    ++scopename;
  }
  
  return scope;
}


PARAMS_HANDLER::SCOPE *
PARAMS_HANDLER::get_curr_set_scope () 
{
  SCOPE *scope = &_set;
  vector<string>::iterator scopename = _scope_nesting.begin();
  
  while (scopename != _scope_nesting.end()) {
    scope = scope->_scopes[*scopename];
    ++scopename;
  }
  
  return scope;
}


const PARAMS_HANDLER::SCOPE *
PARAMS_HANDLER::get_curr_set_scope () const 
{
  PARAMS_HANDLER::SCOPE *scope = const_cast<SCOPE*>(&_set); 
  vector<string>::const_iterator scopename = _scope_nesting.begin();
  
  while (scopename != _scope_nesting.end()) {
    scope = scope->_scopes[*scopename];
    ++scopename;
  }
  
  return scope;
}

PARAMS_HANDLER::SCOPE::SCOPE (string name)
{
  _name = name;
}

PARAMS_HANDLER::SCOPE::~SCOPE ()
{
  _param_map.clear ();
  map<string, SCOPE*>::iterator p;
  for (p=_scopes.begin(); p!=_scopes.end(); ++p)
    delete p->second;
  _scopes.clear ();
}


list <string>
PARAMS_HANDLER::get_scope_names () const
{
  const PARAMS_HANDLER::SCOPE *pc = get_curr_set_scope ();
  
  map<string, SCOPE*>::const_iterator ptr;
  list<string> l;
  ptr = pc->_scopes.begin ();
  while (ptr != pc->_scopes.end()) {
    SCOPE *s = ptr->second;
    l.push_back (s->_name);
    ptr++;
  }
  return l;
}


list <string>
PARAMS_HANDLER::get_param_names () const
{
  const PARAMS_HANDLER::SCOPE *pc = get_curr_set_scope ();
  
  SCOPE::param_map_t::const_iterator ptr;
  ptr = pc->_param_map.begin();
  list<string> l;
  while (ptr != pc->_param_map.end()) {
    l.push_back (ptr->first);
    ptr++;
  }
  return l;
}
