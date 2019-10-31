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
#ifndef ECELL_STEADY_DIFFUSION_H
# define ECELL_STEADY_DIFFUSION_H

#include "ecell.h"
#include "select_les.h"
using STEADY_DIFFUSION::LES;
#include "field.h"
#include "ofactory.h"
#include "mat_diffusion_iso.h"

/**
   The ecells of this class are used by the PROTO_STEADY_DIFFUSION
   protocol.
*/
class ECELL_STEADY_DIFFUSION : public ECELL {

 public: // declarations ////////////////////////////////////////////

  /**
     This is the type of the function that specializations of ECELL_STEADY_DIFFUSION
     need to register with this class so that they can be manufactured
     by the factory of this class.
  */
  typedef ECELL_STEADY_DIFFUSION *(*make_ecell_func) (GCELL *gcell);

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
  static ECELL_STEADY_DIFFUSION *make_ecell (GCELL *gcell, string implementation);

 protected:
  
  /**
     Call this function to register a make function for the ecell
     based on the gcell type specified, and of given implementation.
  */
  static bool register_make_func (make_ecell_func make, string gcell_type, string implementation);

 public: // object functions ////////////////////////////////////////

  /**
   */
  ECELL_STEADY_DIFFUSION (GCELL *gcell) : ECELL(gcell) {
    _mat = 0; _dofmappers.clear ();
  }
  /**
     Assemble the conductivity matrix.
   */
  virtual bool assemble_conductivity_matrix (FIELD_SCALAR *phi, LES *les, FIELD_VECTOR *geometry) = 0;
  /**
     Assemble the source terms.
   */
  virtual bool assemble_source_terms (FIELD_SCALAR *phi, LES *les, FIELD_VECTOR *geometry) = 0;
  /**
     Attach the evaluation cell to the protocol.
     This includes:
     (a) get all needed parameters from the database;
     (b) establish a binding of the degrees of freedom at the nodes
         interacting at the ecell to the LES.
     Arguments: database, les (linear equation solver), field name (what name you wish
     the solver to know the field under), field.
   */
  void attach (class PROTO_STEADY_DIFFUSION *proto, string ggpath);
  /**
     Return the vector of dofmappers.
   */
  vector <DOFMAPPER_SCALAR *> dofmappers () { return _dofmappers; }
  
  /**
   */
  virtual ~ECELL_STEADY_DIFFUSION () {};

  /**
     Set the name of material this cell should use.
  */
  void set_mat (MAT_DIFFUSION_ISO *mat) { _mat = mat; }
  /**
     Get the material the ecell should use.
     */
  MAT_DIFFUSION_ISO *mat () { return _mat; }

 private: // class objects  /////////////////////////////////////////

  static OFACTORY<ECELL_STEADY_DIFFUSION, make_ecell_func, GCELL *> *_ecell_factory;

 private: // object data ////////////////////////////////////////////

  vector <DOFMAPPER_SCALAR *>    _dofmappers;
  MAT_DIFFUSION_ISO             *_mat;
  
};

#endif
