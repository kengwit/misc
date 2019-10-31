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
#ifndef QUANT_EVALUATOR_H
# define QUANT_EVALUATOR_H

#include <list>
#include "algo.h"
#include "gmesh.h"
#include "ecell.h"

/** 
Interface class. Needs to be implemented by a particular protocol.
*/
class QUANT_EVALUATOR {

 public:
  
  static double INVALID_DATA;

 public: // object functions ////////////////////////////////////////
  
  /**
     The quantity evaluator for evaluation cell e, and 
     quantities given by the list of the strings quants.
  */
  QUANT_EVALUATOR (ECELL *e,list<string>quants){
    for(list<string>::iterator i =quants.begin(); i != quants.end(); i++) _quants.push_back(*i);  
    _ecell =e;
  } 
  virtual ~QUANT_EVALUATOR () { _quants.clear(); };
  /**
     Evaluates the location of the point given by the parametric coordinates 
     inside the evaluation cell given to the constructor. 
  */
  virtual FIXED_VECTOR<3> loc (POINT &param_loc) = 0; 
  /**
     Evaluates the quantities at the point given by the parametric coordinates
     inside the evaluation cell given to the constructor. 
  */
  virtual vector <double> quants (POINT &param_loc) = 0;
  list <string> quants () { return _quants;}
  ECELL *ecell() {return _ecell;}
  
  
 private:
  
  list<string> _quants;  
  ECELL       *_ecell;
  
};


#endif
