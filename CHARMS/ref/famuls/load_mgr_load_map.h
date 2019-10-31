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
#ifndef LOAD_MGR_LOAD_MAP_H
# define LOAD_MGR_LOAD_MAP_H

/**
   This class is used to hold instances of loads.
   They may be then retrieved by name.  Protocols would use
   these maps to provide ecells with loaderial handles.
*/
class LOAD_MGR_LOAD_MAP {

 public:

  LOAD_MGR_LOAD_MAP (DB *db) : _db (db) {
    _load_map.clear ();
  }

  /**
     Return the loaderial by name.  It is either found in the map,
     or it is first created and then returned.
  */
  LOAD *load (string loadtype, string loadname) {
    //cerr << "Request to make/find " << loadname << " (" << loadtype << ")" << endl;
    map <string, LOAD *>::iterator i = _load_map.find (loadname);
    if (i == _load_map.end ()) {
      //cerr << "Making " << loadtype << " " << loadname << endl;
      LOAD *m = LOAD::make_load (_db, loadtype, loadname);
      _load_map.insert (map <string, LOAD *>::value_type (loadname, m));
      return m;
    } else {
      //cerr << "Found type " << i->second->type_name() << endl;
      return i->second;
    }
  }

 private:

  DB                       *_db;
  map <string, LOAD *>      _load_map;

};

#endif
