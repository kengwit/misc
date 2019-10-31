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
#include "mat_viscel_iso.h"
extern "C" {
#include "tokensP.h"
}

const string MAT_VISCEL_ISO::TYPE_NAME = MAT_VISCEL::TYPE_NAME + " " +  "viscel_iso";

static MAT *
make (pair <DB *, string> p)
{
  return (new MAT_VISCEL_ISO (p.first, p.second));
}

bool
MAT_VISCEL_ISO::register_make_func ()
{
   tokens_parser _parser = tokens_new_parser ();
   const char *s = MAT_VISCEL_ISO::TYPE_NAME.c_str();
   if (tokens_parse_line (_parser, (char *) s)) {
     for (int tok = 1; tok <= tokens_total_of_tokens (_parser); tok++) {
       //       cerr << "Registering " << string(tokens_token_as_string (_parser, tok)) << endl;
       bool suc = MAT::register_make_func (make, string(tokens_token_as_string (_parser, tok)));
       if (!suc) return false;
     }
   }
   tokens_delete_parser (_parser);
  return true;
  //  return MAT::register_make_func (make, (MAT_ELASTICITY::TYPE_NAME + " " + MAT_ELASTICITY_ISO::TYPE_NAME));
}

MAT_VISCEL_ISO::MAT_VISCEL_ISO (DB *db, string name) : MAT_VISCEL (db, name)
{
  string path = "materials/" + this->name();
  string expr[1];
    if (db->param_defined (path + "/lambda")) {
    expr[0] = db->DB_GET_STRING (path + "/lambda");
  } else {
    if (db->param_defined (path + "/E") &&
        db->param_defined (path + "/nu")) {
      double E  = db->DB_GET_DOUBLE (path + "/E");
      double nu = db->DB_GET_DOUBLE (path + "/nu");
      double lambda = E * nu / (1 + nu) / (1 - 2*nu);
      char buf[512];
      sprintf (buf, "%g", lambda);
      expr[0] = buf;
    } else {
      CHECK (0, EXCEPTION_BAD_VALUE,;);
    }
  }
  _lambda = FUNC<1> (expr);
  if (db->param_defined (path + "/mu")) {
    expr[0] = db->DB_GET_STRING (path + "/mu");
  } else {
    if (db->param_defined (path + "/E") &&
        db->param_defined (path + "/nu")) {
      double E  = db->DB_GET_DOUBLE (path + "/E");
      double nu = db->DB_GET_DOUBLE (path + "/nu");
      double mu = E / 2 / (1 + nu);
      char buf[512];
      sprintf (buf, "%g", mu);
      expr[0] = buf;
    } else {
      CHECK (0, EXCEPTION_BAD_VALUE,;);
    }
  }
  _mu = FUNC<1> (expr);

  if (db->param_defined (path + "/alpha")) {
    expr[0] = db->DB_GET_STRING (path + "/alpha");
    _alpha = FUNC<1> (expr);
  }

  expr[0] = "0";
  if (db->param_defined (path + "/viscosity")) {
    expr[0] = db->DB_GET_STRING (path + "/viscosity");
  }
  _viscosity = FUNC<1> (expr);
  
  expr[0] = "1";
  if (db->param_defined (path + "/rho")) {
    expr[0] = db->DB_GET_STRING (path + "/rho");
  }
  _rho = FUNC<1> (expr);
  
}

bool
MAT_VISCEL_ISO::mat_stiffness (POINT &at, double C_mat[1][1])
{
  double lambda = _lambda (at)(0);
  double mu = _mu (at)(0);
#undef C
#define C(I,J,V) { C_mat[I][J] = V; }
  C(0,0,lambda + 2*mu); 
#undef C
  return true;
}

bool
MAT_VISCEL_ISO::mat_stiffness (POINT &at, double C_mat[3][3])
{
  double lambda = _lambda (at)(0);
  double mu = _mu (at)(0);
#undef C
#define C(I,J,V) { C_mat[I][J] = V; }
  C(0,0,lambda + 2*mu); C(0,1,lambda);        C(0,2,0);  
  C(1,0,lambda);        C(1,1,lambda + 2*mu); C(1,2,0);  
  C(2,0,0);             C(2,1,0);             C(2,2,mu); 
#undef C
  return true;
}

bool
MAT_VISCEL_ISO::mat_stiffness (POINT &at, double C_mat[6][6])
{
  double lambda = _lambda (at)(0);
  double mu = _mu (at)(0);
#undef C
#define C(I,J,V) { C_mat[I][J] = V; }
  C(0,0,lambda + 2*mu); C(0,1,lambda);        C(0,2,lambda);        C(0,3,0);  C(0,4,0);  C(0,5,0);
  C(1,0,lambda);        C(1,1,lambda + 2*mu); C(1,2,lambda);        C(1,3,0);  C(1,4,0);  C(1,5,0);
  C(2,0,lambda);        C(2,1,lambda);        C(2,2,lambda + 2*mu); C(2,3,0);  C(2,4,0);  C(2,5,0);
  C(3,0,0);             C(3,1,0);             C(3,2,0);             C(3,3,mu); C(3,4,0);  C(3,5,0); 
  C(4,0,0);             C(4,1,0);             C(4,2,0);             C(4,3,0);  C(4,4,mu); C(4,5,0); 
  C(5,0,0);             C(5,1,0);             C(5,2,0);             C(5,3,0);  C(5,4,0);  C(5,5,mu); 
#undef C
  return true;
}

