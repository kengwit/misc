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
#ifndef ECELL_SILO_H
# define ECELL_SILO_H

#if defined (SILO) && SILO

#include "ecell.h"
#include "field.h"
#include "silo_writer.h"
#include "evalpt.h"

template <int NUM_COMPONENTS>
class ECELL_SILO : public ECELL {

 public: // declarations ////////////////////////////////////////////

  typedef ECELL_SILO *(*make_ecell_func) (GCELL *gcell);

 public: // class functions  ////////////////////////////////////////

  /**
     Call this function to register a constructor for the ecell
     based on the gcell type specified, and of given implementation.
  */
  static bool register_make_func (make_ecell_func make, string gcell_type, string implementation);
  static ECELL_SILO *make_ecell (GCELL *gcell, string implementation);

 public: // object functions ////////////////////////////////////////

  /**
   */
  ECELL_SILO (GCELL *gcell) : ECELL (gcell) {}
  /**
     Write the SILO record for the ecell to the IO buffer.
     The cell will have the geometry given on input, and the field
     to display on the geometry is also specified.  To evaluate the
     field, the geometry eval_on_geometry is used.
   */
  virtual void stream_silo_data (SILO_WRITER *silo_writer,
                                 FIELD_VECTOR *display_on_geometry,
                                 FIELD_VECTOR *eval_on_geometry,
                                 FIELD<NUM_COMPONENTS> *field);
  /**
     Each implementation has to be able to say which shape
     (in terms of the SILO terminology) it is (hex, tet, polyhedron, ...).
   */
  virtual SILO_WRITER::SILO_UCD_ZONE_SHAPE shape () const = 0;

  /**
   */
  virtual ~ECELL_SILO () {};

 private: // class objects  /////////////////////////////////////////

  static OFACTORY<ECELL_SILO, make_ecell_func, GCELL *> *_ecell_factory;

 private: // object data ////////////////////////////////////////////

};

template <int NUM_COMPONENTS>
OFACTORY<ECELL_SILO<NUM_COMPONENTS>, ECELL_SILO<NUM_COMPONENTS>::make_ecell_func, GCELL *> * ECELL_SILO<NUM_COMPONENTS>::_ecell_factory = 0;

template <int NUM_COMPONENTS>
bool
ECELL_SILO<NUM_COMPONENTS>::register_make_func (ECELL_SILO::make_ecell_func make, string gcell_type, string implementation)
{
  if (_ecell_factory == 0) {
    ECELL_SILO::_ecell_factory = new OFACTORY<ECELL_SILO, ECELL_SILO::make_ecell_func, GCELL* >;
  }
  return (ECELL_SILO::_ecell_factory)->register_make_func (make, gcell_type, implementation);
}

template <int NUM_COMPONENTS>
ECELL_SILO<NUM_COMPONENTS> *
ECELL_SILO<NUM_COMPONENTS>::make_ecell (GCELL *gcell, string implementation)
{
  CHECK (ECELL_SILO<NUM_COMPONENTS>::_ecell_factory, EXCEPTION_NULL_PTR,;);
  return (ECELL_SILO<NUM_COMPONENTS>::_ecell_factory)->make (gcell,
                                                               gcell->type_name (),
                                                               implementation);
}

template <int NUM_COMPONENTS>
void
ECELL_SILO<NUM_COMPONENTS>::stream_silo_data (SILO_WRITER *silo_writer,
                                              FIELD_VECTOR *display_on_geometry,
                                              FIELD_VECTOR *eval_on_geometry,
                                              FIELD<NUM_COMPONENTS> *field) {
}

#endif // SILO enabled?

#endif
