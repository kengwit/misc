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
#ifndef MAT_H
# define MAT_H

#include "ofactory.h"
#include "db.h"
#include "point.h"

class MAT {

 private:

 public: // declarations ////////////////////////////////////////////

  /**
     This is the type of the function that specializations of MAT
     need to register with this class so that they can be manufactured
     by the factory of this class.
  */
  typedef MAT *(*make_mat_func) (pair <DB *, string> arg);

 public: // class functions  ////////////////////////////////////////

  /**
     Clients of this class may call this function to produce instances
     of the derived classes.  They need to pass the database, the material type (string),
     the material name (string).
    */
  static MAT *make_mat (DB *db, string mattype, string matname);

 protected:
  
  /**
     Call this function to register a make function for the material
     based on the name, and of given implementation.
  */
  static bool register_make_func (make_mat_func make, string mattype);

 public:

  MAT (string name) : _name(name) { }
  virtual ~MAT () {}
  
  string name () const { return _name; }
  /**
     This is an abstract class.  The type has to be supplied by the
     specialization.
  */
  virtual string type_name () const = 0;
  virtual double var (string name, POINT &at) = 0;

 private: // class objects  /////////////////////////////////////////

  static OFACTORY<MAT, make_mat_func, pair <DB *, string> > *_mat_factory;

 private:

  string _name;

};

#endif

