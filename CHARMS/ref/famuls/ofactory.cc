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
#include "ofactory.h"
#include <stdio.h>
#include <stdlib.h>

#if 0

class OBJ {

 public:

  OBJ::OBJ () { _gcell_type; _implementation; }
  OBJ::OBJ (string gcell_type, string implementation) { _gcell_type = string (gcell_type); _implementation = string (implementation); }
  void OBJ::display () { cout << "factory object: " << _gcell_type << "; " << _implementation << "\n"; }
  string gcell_type () { return _gcell_type; }
  string implementation () { return _implementation; }

 private:

  string _gcell_type;
  string _implementation;

};

OBJ *make_h20_default (void *)
{
  return new OBJ (string ("h20"), string ("default"));
}

OBJ *make_h8_default (void *)
{
  return new OBJ (string ("h8"), string ("default"));
}

OBJ *make_h8_new (void *)
{
  return new OBJ (string ("h8"), string ("new"));
}

OBJ *make_t3_default (void *)
{
  return new OBJ (string ("t3"), string ("default"));
}

typedef OBJ *(*MAKE_FUNC) (void *ignore);
#define IGNORE 0

int
main (int argc, char **argv)
{
  OFACTORY <OBJ, MAKE_FUNC, void *> ofactory;

  if (!ofactory.register_make_func (&make_h20_default, string ("h20"), string ("default"))) {
    cerr << "Not registered (l." << __LINE__ << ")\n";
    exit (1);
  }
  if (!ofactory.register_make_func (&make_h8_default, string ("h8"), string ("default"))) {
    cerr << "Not registered (l." << __LINE__ << ")\n";
    exit (1);
  }
  if (!ofactory.register_make_func (&make_h8_new, string ("h8"), string ("new"))) {
    cerr << "Not registered (l." << __LINE__ << ")\n";
    exit (1);
  }
  if (!ofactory.register_make_func (&make_t3_default, string ("t3"), string ("default"))) {
    cerr << "Not registered (l." << __LINE__ << ")\n";
    exit (1);
  }
  

  {
    OBJ *o = ofactory.make (IGNORE, string ("h20"), string ("default"));
    if (!o) cerr << "Not found (l." << __LINE__ << ")\n";
    else     {
      o->display ();
      delete o;
    }
  }

  {
    OBJ *o = ofactory.make (IGNORE, string ("h20"), string ("total"));
    if (!o) cerr << "Not found (l." << __LINE__ << ")\n";
    else     {
      o->display ();
      delete o;
    }
  }
  
  {
    OBJ *o = ofactory.make (IGNORE, string ("h8"), string ("new"));
    if (!o) cerr << "Not found (l." << __LINE__ << ")\n";
    else     {
      o->display ();
      delete o;
    }
  }
  
  {
    OBJ *o = ofactory.make (IGNORE, string ("h8"), string ("eulhex"));
    if (!o) cerr << "Not found (l." << __LINE__ << ")\n";
    else     {
      o->display ();
      delete o;
    }
  }
  
  {
    OBJ *o = ofactory.make (IGNORE, string ("h8"), string ("default"));
    if (!o) cerr << "Not found (l." << __LINE__ << ")\n";
    else     {
      o->display ();
      delete o;
    }
  }
  
  {
    OBJ *o = ofactory.make (IGNORE, string ("funny"), string ("doesn't exist"));
    if (!o) cerr << "Not found (l." << __LINE__ << ")\n";
    else     {
      o->display ();
      delete o;
    }
  }
  
  {
    OBJ *o = ofactory.make (IGNORE, string ("t3"), string ("default"));
    if (!o) cerr << "Not found (l." << __LINE__ << ")\n";
    else     {
      o->display ();
      delete o;
    }
  }
    
  {
    OBJ *o = ofactory.make (IGNORE, string ("t3"), string ("default"));
    if (!o) cerr << "Not found (l." << __LINE__ << ")\n";
    else     {
      o->display ();
      delete o;
    }
  }

  {
    OBJ *o = ofactory.make (IGNORE, string ("t3"), string ("default"));
    if (!o) cerr << "Not found (l." << __LINE__ << ")\n";
    else     {
      o->display ();
      delete o;
    }
  }

  exit (0);
}

#endif
