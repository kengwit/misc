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
#ifndef SMART_HANDLE_H
#  define SMART_HANDLE_H

#include "ref_counter.h"

/**
   Smart pointer (handle) class.  The pointers are reference-counted.
   When the last pointer goes out of scope, it will delete the referenced
   (pointed-to) object. The smart pointer will point either to a valid object, 
   or it will be null.
   
   NOTE: the object pointed to must be allocated on the heap!
*/

template <typename TYP>
class SMART_HANDLE {
  
 public: // object methods ///////////////////////////////////////
  
  /**
     Default constructor.
  */
  SMART_HANDLE (TYP *p = NULL) : _ptr(p), _ref_counter(NULL) {
    if (_ptr != NULL) {
      _ref_counter = new REF_COUNTER();
    }
  }
  
  /**
     Template copy constructor from OTYP that can be implicitly cast to type TYP
  */
  template<class OTYP>
    SMART_HANDLE (const SMART_HANDLE<OTYP>& smh) 
    : _ptr (smh.get_ptr()) ,
    _ref_counter (smh.get_ref_counter()) 
    {
      adj_ref_counter();
    }
  
  /**
   */
  SMART_HANDLE (const SMART_HANDLE<TYP>& smh) 
    : _ptr (smh.get_ptr()) , _ref_counter (smh.get_ref_counter()) 
    {
      adj_ref_counter();
    }

  /**
     Destructor.  Deletes a reference, which could mean
     the object is no longer referenced by anybody, and will
     get deleted.
  */
  ~SMART_HANDLE (void) { reset_ptr(); }

  /**
     Template assignment operator that takes type OTYP which can be
     implicitly cast to type TYP
  */
  template<class OTYP>
    SMART_HANDLE& operator= (const SMART_HANDLE<OTYP>& rhs) 
    {
      retarget (rhs.get_ptr(), rhs.get_ref_counter()) ;
      return (*this) ;
    }
  
  /**
     Equal? Type OTYP compared with a compatible type TYP.
  */
  template<class OTYP>
    inline bool operator== (const SMART_HANDLE<OTYP>& rhs)  const
    {
      return (_ptr == rhs.get_ptr()) ;
    }
  
  /**
     Compares a type OTYP to a compatible type TYP
  */
  template<class OTYP>
    inline bool operator!= (const SMART_HANDLE<OTYP>&rhs)  const
    {
      return (_ptr != rhs.get_ptr()) ;
    }
  
  /**
     Assignment that resets the pointer value.
  */
  SMART_HANDLE& operator= (const SMART_HANDLE<TYP>& rhs) {
    retarget (rhs.get_ptr(), rhs.get_ref_counter()) ;
    return (*this) ;
  }

  /**
     Assignment.
  */
  SMART_HANDLE& operator= (TYP* p) {
    reset_ptr();
    _ptr = p;
    if (_ptr != NULL) {
      _ref_counter = new REF_COUNTER();
    } else {
      _ref_counter = NULL;
    }
    return (*this) ;
  }

  /**
     Allow access to the data pointer.
  */
  inline TYP& operator* (void)  const throw (EXCEPTION_NULL_PTR) {
    CHECK (_ptr != NULL, EXCEPTION_NULL_PTR, ;);
    return (*_ptr) ;
  }

  /**
     Allow access to the data pointer.
  */
  inline TYP* operator-> (void)  const throw (EXCEPTION_NULL_PTR) {
    if (_ptr == NULL) 
      cout<<"\n Pointer is NULL \n";  
    CHECK (_ptr != NULL, EXCEPTION_NULL_PTR,;);
    return _ptr;
  }

  /**
     Comparing with 0: Is the pointer NULL?
  */
  inline bool operator== (int zero)  const throw (EXCEPTION_BAD_VALUE) {
    CHECK (zero == 0, EXCEPTION_BAD_VALUE,;);
    return (_ptr == NULL) ;
  }

  /**
     Comparing with a non-null pointer.
  */
  inline bool operator!= (int zero)  const throw (EXCEPTION_BAD_VALUE) {
    CHECK (zero == 0, EXCEPTION_BAD_VALUE,;);
    return (_ptr != NULL) ;
  }
  
  /**
     Reset the pointer.  One reference is removed,
     which could cause the data pointer to get deleted.
  */
  inline void reset_ptr (void) {
    if (_ref_counter != NULL) {
      bool del = _ref_counter->remove_ref();
      if (del) {
        delete _ptr;
        delete _ref_counter;
      }
      _ref_counter = NULL;
      _ptr = NULL;
    }
  }

  inline TYP* get_ptr (void)  const {
    return _ptr;
  }
  inline REF_COUNTER* get_ref_counter (void)  const {
    return _ref_counter;    
  }

  /** Forget the pointer.
      The pointer must be the one currently linked to the handle:
      otherwise nothing happens.
   */
  void forget_me (TYP *ptr) {
    if (ptr == _ptr)
      if (_ref_counter != NULL) {
        bool del = _ref_counter->remove_ref();
        if (del) { // Note: the ptr is *not* deleted!
          delete _ref_counter;
        }
        _ref_counter = NULL;
        _ptr = NULL;
      }
  }
  
 private: // object data
  
  TYP*         _ptr;
  REF_COUNTER* _ref_counter;
  
 private: // object methods ///////////////////////////////////

  inline void adj_ref_counter (void) {
    if (_ref_counter != NULL) {
      if (!_ref_counter->add_ref()) {
        _ref_counter = NULL;
        _ptr = NULL;
      }
    }
  }
  
  inline void retarget (TYP* p, REF_COUNTER* cntr) {
    if (_ref_counter != cntr)  {
      reset_ptr();
      _ptr = p;
      _ref_counter = cntr;
      adj_ref_counter();
    }
  }
  
};

#endif
