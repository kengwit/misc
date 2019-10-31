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
#ifndef MAT_MGR_MAT_MAP_H
# define MAT_MGR_MAT_MAP_H

/**
   This class is used to hold instances of materials.
   They may be then retrieved by name.  Protocols would use
   these maps to provide ecells with material handles.
*/
class MAT_MGR_MAT_MAP {

 public:

  MAT_MGR_MAT_MAP (DB *db) : _db (db) {
    _mat_map.clear ();
  }

  /**
     Return the material by name.  It is either found in the map,
     or it is first created and then returned.
  */
  MAT *mat (string mattype, string matname) {
    //cerr << "Request to make/find " << matname << " (" << mattype << ")" << endl;
    map <string, MAT *>::iterator i = _mat_map.find (matname);
    if (i == _mat_map.end ()) {
      //cerr << "Making " << mattype << " " << matname << endl;
      MAT *m = MAT::make_mat (_db, mattype, matname);
      _mat_map.insert (map <string, MAT *>::value_type (matname, m));
      return m;
    } else {
      //cerr << "Found type " << i->second->type_name() << endl;
      return i->second;
    }
  }

 private:

  DB                       *_db;
  map <string, MAT *>       _mat_map;

};

#endif
