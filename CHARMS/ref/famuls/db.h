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
#ifndef DB_H
# define DB_H

#include "params.h"

/**
   This class implements a very simple ASCII parameter file database.
   Parameters are accessed as if they were parts of a class/object hierarchy.

   Example (`!' starts a comment):
   
<pre>
   obj algorithms               ! class name
    obj steady_diffusion        ! class name
      obj adsq                  ! an instance
        param gmesh = m1        ! parameter: mesh name
        param gsubmeshes = sm1  ! parameter: work with submeshes of the above mesh
        obj gcell_groups        ! parameter: 
          obj gg1               ! specification of group gg1
            param implementation = default ! parameter
            param material = generic_heat_diffusion ! parameter
          end
        end
        param les = default_les ! parameter: linear equation solver
        obj essential_bcs       ! specification of EBC's
          param ebc_file = adsq.ebc ! parameter: EBC file
        end
        obj adapt               ! specification of adaptivity
          param adaptive = true
        end
        param print_sol = false ! parameter: printouts?
      end
    end
</pre>
*/
class DB {

 public:

  DB ();
  ~DB ();

  bool load_file (string file);
  bool load_string (string s);
  bool param_defined (const string param_path);
  string get_string_from_path (char *file, int line, const string param_path);
  double get_double_from_path (char *file, int line, const string param_path);
  long int get_integer_from_path (char *file, int line, const string param_path);
  bool get_bool_from_path (char *file, int line, const string param_path);
  list <string> get_string_list_from_path (char *file, int line, const string param_name);
  list <long int> get_integer_list_from_path (char *file, int line, const string param_name);
  list <double> get_double_list_from_path (char *file, int line, const string param_name);

 private:
  
  PARAMS_HANDLER *_params_handler;
  
};

#define DB_GET_STRING(param_path) get_string_from_path (__FILE__, __LINE__, param_path)
#define DB_GET_DOUBLE(param_path) get_double_from_path (__FILE__, __LINE__, param_path)
#define DB_GET_INTEGER(param_path) get_integer_from_path (__FILE__, __LINE__, param_path)
#define DB_GET_BOOL(param_path) get_bool_from_path (__FILE__, __LINE__, param_path)
#define DB_GET_STRING_LIST(param_path) get_string_list_from_path (__FILE__, __LINE__, param_path)
#define DB_GET_DOUBLE_LIST(param_path) get_double_list_from_path (__FILE__, __LINE__, param_path)
#define DB_GET_INTEGER_LIST(param_path) get_integer_list_from_path (__FILE__, __LINE__, param_path)


#endif
