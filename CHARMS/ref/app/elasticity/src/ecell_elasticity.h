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
#ifndef ECELL_ELASTICITY_H
# define ECELL_ELASTICITY_H

#include "ecell.h"
#include "select_les.h"
using ELASTICITY::LES;
#include "field.h"
#include "ofactory.h"
#include "func.h"
#include "mat_elasticity.h"
#include "load_on_gcell_body.h"
#include "load_on_gcell_traction.h"

/**
   The ecells of this class are used by the PROTO_ELASTICITY
   protocol.
*/
class ECELL_ELASTICITY : public ECELL {

 public: // declarations ////////////////////////////////////////////

  /**
     This is the type of the function that specializations of ECELL_ELASTICITY
     need to register with this class so that they can be manufactured
     by the factory of this class.
  */
  typedef ECELL_ELASTICITY *(*make_ecell_func) (GCELL *gcell);

 public: // class functions  ////////////////////////////////////////

  /**
     Clients of this class may call this function to produce instances
     of the derived classes.  They need to pass an instance of a gcell, and
     a string indicating which implementation of the ecell should be created.
     In other words, given (a) the type of the gcell, and (b) a descriptive name
     (string) that distinguishes among different implementations of the same behaviour,
     the factory of this class will produce an instance of this class (actually,
     a descendant of this class that overrides the required methods).
   */
  static ECELL_ELASTICITY *make_ecell (GCELL *gcell, string implementation);

 protected:
  
  /**
     Call this function to register a make function for the ecell
     based on the gcell type specified, and of given implementation.
  */
  static bool register_make_func (make_ecell_func make, string gcell_type, string implementation);

 public: // object functions ////////////////////////////////////////

  /**
   */
  ECELL_ELASTICITY (GCELL *gcell) : ECELL(gcell) {
    _mat = 0; _dofmappers.clear ();
    _body_load = 0;
    _traction_load = 0;
  }
  /**
     Assemble the stiffness matrix.
   */
  virtual bool assemble_stiffness_matrix (FIELD_VECTOR *u, LES *les, FIELD_VECTOR *geometry) = 0;
  /**
     Assemble the load due to prescribed displacements.
   */
  virtual bool assemble_prescribed_u_load (FIELD_VECTOR *u, LES *les, FIELD_VECTOR *geometry) = 0;
  /**
     Assemble the load due to prescribed body loads.
   */
  virtual bool assemble_body_load (FIELD_VECTOR *u, LES *les, FIELD_VECTOR *geometry) = 0;
  /**
     Assemble the load due to prescribed tractions.
   */
  virtual bool assemble_tractions_load (FIELD_VECTOR *u, LES *les, FIELD_VECTOR *geometry) = 0;
  /**
     Attach the evaluation cell to the protocol.
     This includes:
     (a) get all needed parameters from the database;
     (b) establish a binding of the degrees of freedom at the nodes
         interacting at the ecell to the LES.
     Arguments: database, les (linear equation solver), field name (what name you wish
     the solver to know the field under), field.
   */
  void attach (class PROTO_ELASTICITY *proto, string ggpath);
  /**
     Return the vector of dofmappers.
   */
  vector <DOFMAPPER_VECTOR *> dofmappers () { return _dofmappers; }
  
  /**
   */
  virtual ~ECELL_ELASTICITY () {};

  /**
     Set the name of material this cell should use.
  */
  void set_mat (MAT_ELASTICITY *mat) { _mat = mat; }
  /**
     Get the material the ecell should use.
     */
  MAT_ELASTICITY * mat () const { return _mat; }
  /**
     Get the body load (could be null, in which case no body load is defined).
  */
  LOAD_ON_GCELL_BODY *body_load () const { return _body_load; }
  /**
     Get the traction load (could be null, in which case no traction load is defined).
  */
  LOAD_ON_GCELL_TRACTION *traction_load () const { return _traction_load; }

 private: // class objects  /////////////////////////////////////////

  static OFACTORY<ECELL_ELASTICITY, make_ecell_func, GCELL *> *_ecell_factory;

 private: // object data ////////////////////////////////////////////

  vector <DOFMAPPER_VECTOR *>    _dofmappers;
  MAT_ELASTICITY                *_mat;
  LOAD_ON_GCELL_BODY            *_body_load;
  LOAD_ON_GCELL_TRACTION        *_traction_load;
  
};

#endif
