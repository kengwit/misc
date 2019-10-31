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
#ifndef LOAD_H
# define LOAD_H

#include "ofactory.h"
#include "db.h"
#include "point.h"

class LOAD {

 private:

 public: // declarations ////////////////////////////////////////////

  /**
     This is the type of the function that specializations of LOAD
     need to register with this class so that they can be manufactured
     by the factory of this class.
  */
  typedef LOAD *(*make_load_func) (pair <DB *, string> arg);

 public: // class functions  ////////////////////////////////////////

  /**
     Clients of this class may call this function to produce instances
     of the derived classes.  They need to pass the database, the loaderial type (string),
     the loaderial name (string).
    */
  static LOAD *make_load (DB *db, string loadtype, string loadname);

 protected:
  
  /**
     Call this function to register a make function for the loaderial
     based on the name, and of given implementation.
  */
  static bool register_make_func (make_load_func make, string loadtype);

 public:

  LOAD (string name) : _name(name) { }
  virtual ~LOAD () {}
  
  string name () const { return _name; }
  /**
     This is an abstract class.  The type has to be supplied by the
     specialization.
  */
  virtual string type_name () const = 0;
  virtual double var (string name, POINT &at) = 0;
  virtual double var (string name, POINT4 &at) = 0;

 private: // class objects  /////////////////////////////////////////

  static OFACTORY<LOAD, make_load_func, pair <DB *, string> > *_load_factory;

 private:

  string _name;

};

#endif

