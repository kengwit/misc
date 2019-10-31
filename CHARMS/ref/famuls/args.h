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
#ifndef ARGS_H
# define ARGS_H

#include <string>
extern "C" {
#include "getoptP.h"
}

class ARGS {

 public: // class objects ///////////////////////////////////////////

  static int argc;
  static char **argv;

  static void start_process_options (int argc, char **argv, std::string usagestring);
  static void finish_process_options ();
  static bool get_opt_bool (std::string optstring, bool default_value, std::string helpstring);
  static int get_opt_int (std::string optstring, int default_value, std::string helpstring);
  static double get_opt_double (std::string optstring, double default_value, std::string helpstring);
  static std::string get_opt_string (std::string optstring, std::string default_value, std::string helpstring);

 private: // class data /////////////////////////////////////////////

  static bool        _report_opts;
  static bool        _help;
  static std::string _usagestring;

  static void usage ();
  
};

#endif
