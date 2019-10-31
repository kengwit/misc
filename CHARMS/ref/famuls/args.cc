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
#include "args.h"
#undef min // conflict with std_limits.h
#undef max // conflict with std_limits.h
#include <iostream>

using namespace std;

int ARGS::argc = 0;
char ** ARGS::argv = 0;
bool ARGS::_report_opts = false;
bool ARGS::_help = false;
string ARGS::_usagestring = "";

void
ARGS::usage ()
{
  char *option, *usage;

  cerr << "Usage: " << ARGS::argv[0] << " " << ARGS::_usagestring << " + options below" << endl;
  
  if (ckit_opt_first_recognized (&option, &usage)) {
    cerr << (option != NULL ? option : "?") << " " << (usage  != NULL ? usage  : "?") << endl;
  }
  while (ckit_opt_next_recognized (&option, &usage)) {
    cerr << (option != NULL ? option : "?") << " " << (usage  != NULL ? usage  : "?") << endl;
  }
  exit (0);
}

void 
ARGS::start_process_options (int argc, char **argv, string usagestring)
{
  ARGS::argc = argc;
  ARGS::argv = argv;
  ARGS::_usagestring = usagestring;
  
  ARGS::_report_opts = get_opt_bool ("-report_opts", false, "prints options as they are being set");
  ARGS::_help        = get_opt_bool ("--help", false, "prints usage for each option recognized");
}

bool
ARGS::get_opt_bool (string optstring, bool default_value, string helpstring)
{
  bool result = default_value;
  if (ckit_opt_is_set (ARGS::argc, ARGS::argv,
                       (char *) optstring.c_str(), (char *) helpstring.c_str(), 1)) {
    result = true;
  } 
  if (_report_opts) { cerr << "Set: " << optstring << " to true" << endl; }
  return result;
}

int
ARGS::get_opt_int (string optstring, int default_value, string helpstring)
{
  int result = default_value;
  int how_many_occurrences = ckit_opt_is_set (ARGS::argc, ARGS::argv,
                                              (char *) optstring.c_str(),
                                              (char *) helpstring.c_str(), 1);
  if (how_many_occurrences == 1) {
    if (ckit_opt_has_arg ()) {
      if (!ckit_opt_int_val (&result)) {
        cerr << "Bad argument for " << optstring << endl;
        abort ();
      }
    } else {
      cerr << "Expected argument for " << optstring << endl;
      abort ();
    }
  } else if (how_many_occurrences > 1) {
    cerr << "Only one option per command line is allowed for " << optstring << endl;
    abort ();
  } 
  if (_report_opts) { cerr << "Set: " << optstring << " to " << result << endl; }
  return result;
}

double
ARGS::get_opt_double (string optstring, double default_value, string helpstring)
{
  double result = default_value;
  int how_many_occurrences = ckit_opt_is_set (ARGS::argc, ARGS::argv,
                                              (char *) optstring.c_str(),
                                              (char *) helpstring.c_str(), 1);
  if (how_many_occurrences == 1) {
    if (ckit_opt_has_arg ()) {
      if (!ckit_opt_double_val (&result)) {
        cerr << "Bad argument for " << optstring << endl;
        abort ();
      }
    } else {
      cerr << "Expected argument for " << optstring << endl;
      abort ();
    }
  } else if (how_many_occurrences > 1) {
    cerr << "Only one option per command line is allowed for " << optstring << endl;
    abort ();
  }
  if (_report_opts) { cerr << "Set: " << optstring << " to " << result << endl; }
  return result;
}

string
ARGS::get_opt_string (string optstring, string default_value, string helpstring)
{
  string result = default_value;
  int how_many_occurrences = ckit_opt_is_set (ARGS::argc, ARGS::argv,
                                              (char *) optstring.c_str(),
                                              (char *) helpstring.c_str(), 1);
  if (how_many_occurrences == 1) {
    if (ckit_opt_has_arg ()) {
      char *s = 0;
      if (!ckit_opt_char_val (&s)) {
        cerr << "Bad argument for " << optstring << endl;
        abort ();
      }
      result = s;
    } else {
      cerr << "Expected argument for " << optstring << endl;
      abort ();
    }
  } else if (how_many_occurrences > 1) {
    cerr << "Only one option per command line is allowed for " << optstring << endl;
    abort ();
  } 
  if (_report_opts) { cerr << "Set: " << optstring << " to " << result << endl; }
  return result;
}

void
ARGS::finish_process_options ()
{
  if (ARGS::_help) {
    ARGS::usage ();
  }
}
