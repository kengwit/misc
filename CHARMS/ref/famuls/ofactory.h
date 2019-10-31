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
#ifndef OFACTORY_H
# define OFACTORY_H

#include <map>
#include <string>

template <typename OBJ_TO_MAKE, typename MAKE_FUNC, typename MAKE_ARG>
class OFACTORY {
  
 public:

  OFACTORY () { }
  
  bool register_make_func (MAKE_FUNC make_func, std::string type, std::string implementation) {
    MAKE_FUNC ac = find_make_func (type, implementation);
    if (ac) return false;
    else    return add_make_func (make_func, type, implementation);
  }

  OBJ_TO_MAKE *make (MAKE_ARG arg, std::string type, std::string implementation) {
    MAKE_FUNC c = find_make_func (type, implementation);
    if (!c) return 0;
    else    return (*c) (arg);
  }

 private:

  typedef std::map <std::string, MAKE_FUNC>  implementation_lookup_map;
  typedef std::map <std::string, implementation_lookup_map *>  type_lookup_map;
  
  type_lookup_map _type_lookup_map;
  
  MAKE_FUNC find_make_func (std::string type, std::string implementation) {
    implementation_lookup_map *ilmap = find_implementation_map (type);
    if (ilmap != 0) {
      typename implementation_lookup_map::iterator ilmapi = ilmap->find (implementation);
      if (ilmapi != ilmap->end ()) {
        MAKE_FUNC c = ilmapi->second;
        return c;
      } else
        return 0;
    } else
      return 0;
  }

  implementation_lookup_map *find_implementation_map (std::string type) {
    typename type_lookup_map::iterator gtlmapi = _type_lookup_map.find (type);
    if (gtlmapi != _type_lookup_map.end ()) {
      implementation_lookup_map *ilmap = gtlmapi->second;
      return ilmap;
    } else
      return 0;
  }

  bool add_make_func (MAKE_FUNC make_func, std::string type, std::string implementation) {
    implementation_lookup_map *ilmap = find_implementation_map (type);
    if (ilmap == 0) {
      ilmap = new implementation_lookup_map ();
      if (ilmap == 0) return false;
      (_type_lookup_map).insert (typename type_lookup_map::value_type (type, ilmap));
    }
    (*ilmap).insert (typename implementation_lookup_map::value_type (implementation, make_func));
    return (find_make_func (type, implementation) != 0);
  }
  
};
  
#endif
