
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
#ifndef ENUMERATOR_H
#  define ENUMERATOR_H

template <typename CONTAINER, typename OBJ>
class ENUMERATOR {
  
 public: // object functions ////////////////////////////////////////
  
  /**
    Construct an enumerator for the STL container on input. 
   */
  ENUMERATOR (CONTAINER c) {
    _container = c;
    reset ();
  }
  virtual ~ENUMERATOR (){};
  /**
    Reset the enumerator so that the next element obtained
    from the next() call is the first object in the container.
     */
  void reset () { _current = _container.begin(); }
  /**
    Get the next object from the enumerator.  Zero (0) is returned
    when the enumerator arrives at the end of the associated container.
     */
  virtual OBJ *next () {
    if (_current == _container.end()) return 0;
    else {
      OBJ *result = (*_current);
      _current++;
      return result;
    }
  }
  
 protected:
  
  CONTAINER                          _container;
  typename CONTAINER::const_iterator _current;
  
};

#endif
