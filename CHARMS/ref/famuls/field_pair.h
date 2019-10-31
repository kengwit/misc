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
#ifndef FIELD_PAIR_H
# define FIELD_PAIR_H

#include "algorithm"
#include "list"
#include "bfun.h"
#include "fen.h"
#include "fixed_vector.h"

template <int NUM_COMPONENTS>
class FIELD_PAIR {

 public: // object functions ////////////////////////////////////////

  /**
     Create a field pair (bfun + dofparam).
  */
  FIELD_PAIR (BFUN *bfun, FIXED_VECTOR<NUM_COMPONENTS> initial_value) {
    _bfun        = bfun;
    _dofparam    = initial_value;
    _constrained = 0;
  }

  ~FIELD_PAIR () { if (_constrained != 0) delete [] _constrained; }
  
  /**
     Return the basis function.
  */
  BFUN *bfun () const { return _bfun; }

  /**
     How many components in a dofparam?
  */
  int num_components () const { return NUM_COMPONENTS; }

  /**
     Set the value of the parameter of the field pair.
     If any component of the parameter is constrained (ie. not free),
     such component is unchanged.
  */
  void set_dofparam (FIXED_VECTOR<NUM_COMPONENTS> value) {
    for (int j = 0; j < NUM_COMPONENTS; j++) set_dofparam_comp (j, value[j]);
  }
  /**
     Set the value of the parameter of the field pair.
     If any component of the parameter is unconstrained (ie. free),
     such component is unchanged.
  */
  void set_dofparam_constrained_only (FIXED_VECTOR<NUM_COMPONENTS> value) {
    for (int j = 0; j < NUM_COMPONENTS; j++) 
      set_dofparam_comp_constrained_only (j, value[j]);
  }
  /**
     Set the value of the parameter of the field pair.
     All components are set, irrespectively of whether they are constrained
     or not.
  */
  void set_dofparam_all (FIXED_VECTOR<NUM_COMPONENTS> value) {
    for (int j = 0; j < NUM_COMPONENTS; j++) {
      //CHECK(!isnan(value[j]), EXCEPTION_BAD_VALUE,;); // RAW: expensive, use only for debug
      _dofparam[j] = value[j];
    }
  }

  /**
     Set the value of a single component of the parameter of the field pair.
     If the component is constrained, its value will not change.
  */
  void set_dofparam_comp (int comp, double value) {
     if (_constrained) {
      if (! _constrained[comp])
        _dofparam[comp] = value;
    } else
      _dofparam[comp] = value;
//   cout << "setting[" << comp << "] " << _dofparam[comp] << endl; //RAW:ashi
  }
  /**
     Set the value of a single component of the parameter of the field pair.
     If the component is unconstrained (free), its value will *not* change.
  */
  void set_dofparam_comp_constrained_only (int comp, double value) {
    //CHECK(!isnan(value), EXCEPTION_BAD_VALUE,;); // RAW: expensive, use only for debug
    if (_constrained) {
      if (_constrained[comp])
        _dofparam[comp] = value;
    } 
  }

  /**
     Get the value of the parameter.
  */
  FIXED_VECTOR<NUM_COMPONENTS> dofparam () const {
    // cerr << "getting dofparam of node " << _bfun->fen()->id() << " " << _dofparam[0] << endl;
    return _dofparam;
  }

  /**
     Get the value of a single component of the parameter.
  */
  double dofparam_comp (int comp) const {
    // cerr << "getting[" << comp << "] " << _dofparam[comp] << endl;
    return _dofparam[comp];
  }

  /**
     Any component of the parameter constrained?  If any component is
     constrained, the function returns true; otherwise the function returns
     false.
  */
  bool is_constrained () const {
    if (_constrained == 0) {
      return false;
    } else {
      for (sizet j = 0; j < NUM_COMPONENTS; j++) {
        if (_constrained[j]) return true;
      }
      return false;
    }
  }
  /**
     If any component of the dofparam is constrained,  the corresponding
     constrained[j] element is true; otherwise the values in the array constrained are all false. 
  */
  void constrained (bool constrained[NUM_COMPONENTS]) const {
    if (_constrained == 0) {
      for (int j = 0; j < NUM_COMPONENTS; j++) { constrained[j] = false; }
    } else {
      for (int j = 0; j < NUM_COMPONENTS; j++) {
        constrained[j] = _constrained[j]; 
      }
    }
  }
  /**
     Set value and constraint of a component.
  */
  void set_constrained_comp (int J, bool constrained, double value) {
    _dofparam[J] = value; // set the value
    if (constrained) { // set the constraint
      if (_constrained == 0) { // initialize
        _constrained = new bool[NUM_COMPONENTS];
        for (int j = 0; j < NUM_COMPONENTS; j++) { _constrained[j] = false; }
      }
      _constrained[J] = true; // set the constraint
    } else { // free the component
      if (_constrained) { // any constraints had been set?
        _constrained[J] = false; // release
        bool any_constrained = false; // now check if we have to keep _constrained around
        for (int j = 0; j < NUM_COMPONENTS; j++) {
          if (_constrained[j]) {
            any_constrained = true;
            break;
          }
        }
        if (!any_constrained) { // no constraints are set
          delete [] _constrained;
          _constrained = 0;
        }
      }
    }
  }

 private:

  BFUN                        *_bfun;
  FIXED_VECTOR<NUM_COMPONENTS> _dofparam;
  bool                        *_constrained;

};

#endif
