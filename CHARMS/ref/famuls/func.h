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
#ifndef FUNC_H
# define FUNC_H

#include <string>
#include "func_base.h"
#include "expevalP.h"
#include "tokensP.h"
#include "fixed_vector.h"
#include "point.h"

template <int NUM_COMPONENTS>
class FUNC: public FUNC_BASE {

 public:

  FUNC () : _is_constant(true), _constant_value(0.0) {
    for (sizet j = 0; j < NUM_COMPONENTS; j++) {
      _expr[j] = "";
    }
  }
  
  FUNC (double v) : _is_constant(true), _constant_value(v) {
    for (sizet j = 0; j < NUM_COMPONENTS; j++) {
      _expr[j] = "";
    }
  }
  
  FUNC (string expr[NUM_COMPONENTS]) {
    tokens_parser parser = tokens_new_parser ();
    CHECK (parser, EXCEPTION_NULL_PTR,;);
    _is_constant = true;
    for (sizet j = 0; j < NUM_COMPONENTS; j++) {
      _expr[j] = expr[j];
      if (_expr[j] == "\"\"") _expr[j] = ""; // treat empty strings
      if (_expr[j] != "") { // Try to parse it
        CHECK (tokens_parse_line (parser, (char *) _expr[j].c_str ()), EXCEPTION_BAD_VALUE,;);
        if ((tokens_total_of_tokens (parser) == 1) &&
            (tokens_type_of_token (parser, 1) == TOKENS_NUMBER)) { // Is the token a number?
          _is_constant =  _is_constant &&  true; 
          _constant_value(j) = tokens_token_as_double (parser, 1);
          if ((tokens_total_of_tokens (parser) == 1) &&
              (tokens_token_and_keyword_equiv(tokens_token_as_string(parser,1),"" ))) {
            _expr[j] = "";
          }
        } else if ((tokens_total_of_tokens (parser) == 1) &&
                   (tokens_type_of_token (parser, 1) == TOKENS_STRING)) {
          _is_constant = false;
          _expr[j] = tokens_token_as_string (parser, 1);
          _constant_value(j) = 0.0;
        } else { // Expression: cannot assume constantness
          _is_constant = false;
          _constant_value(j) = 0.0;
        }
      } else {
        _constant_value(j) = 0.0;
      }
    }
    tokens_delete_parser (parser);
  }

  virtual ~FUNC () {
  }

  /**
     Evaluate the function at a point.  If a component of the vector
     function is not defined, zero (0) is returned (and should be ignored,
     of course).
  */
  virtual FIXED_VECTOR<NUM_COMPONENTS> operator() (POINT& at) {
    if (_is_constant) {
      return _constant_value;
    } else {
      FIXED_VECTOR<NUM_COMPONENTS> v(0.0);
      expeval_set_value ("x", at(0));
      expeval_set_value ("y", at(1));
      expeval_set_value ("z", at(2));
      for (sizet j = 0; j < NUM_COMPONENTS; j++) {
        double val = 0;
        if (_expr[j] != "" ) {
          if (expeval_eval ((char *) _expr[j].c_str (), &val) != EXPEVAL_E_OK) {
            cerr << "Expression: " <<_expr[j] << endl;
            CHECK (0, EXCEPTION_BAD_VALUE,;);
          //cerr << "evaluated func " << _expr[j] << " as " << val << endl;
          }
        } 
        v(j) = val;
      }
      return v;
    }
  }

  /**
     Evaluate the function at a point.  If a component of the vector
     function is not defined, zero (0) is returned (and should be ignored,
     of course).
     Here the point is FIXED_VECTOR<4>, which incorporates time in at(3).
  */
  virtual FIXED_VECTOR<NUM_COMPONENTS> operator() (POINT4& at) {
    if (_is_constant) {
      return _constant_value;
    } else {
      FIXED_VECTOR<NUM_COMPONENTS> v(0.0);
      expeval_set_value ("x", at(0));
      expeval_set_value ("y", at(1));
      expeval_set_value ("z", at(2));
      expeval_set_value ("t", at(3));
      
      for (sizet j = 0; j < NUM_COMPONENTS; j++) {
        double val = 0;
        if ( _expr[j] != "") {
          if (expeval_eval ((char *) _expr[j].c_str (), &val) != EXPEVAL_E_OK) {
            cerr << "Expression: " <<_expr[j] << endl;
            CHECK (0, EXCEPTION_BAD_VALUE,;);
           }
        } 
        v(j) = val;
      }
      return v;
    }
  }


  /**
     Is this component of the function given?
  */
  virtual bool given (int j) {
    return (_expr[j] != "");
  }

 private:

  bool                         _is_constant;
  FIXED_VECTOR<NUM_COMPONENTS> _constant_value;
  string                       _expr[NUM_COMPONENTS];

};

typedef FUNC<1> FUNC_SCALAR;
typedef FUNC<3> FUNC_VECTOR;

#endif
