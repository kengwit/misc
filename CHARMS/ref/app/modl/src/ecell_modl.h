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
#ifndef ECELL_MODL_H
# define ECELL_MODL_H

#include "ecell.h"
#include "evs_ooofs.h"
#include "field.h"
#include "ofactory.h"
#include "func.h"
#include "mat_elasticity.h"

/**
   The ecells of this class are used by the PROTO_MODL
   protocol.
*/
class ECELL_MODL : public ECELL {

 public: // declarations ////////////////////////////////////////////

  /**
     This is the type of the function that specializations of ECELL_MODL
     need to register with this class so that they can be manufactured
     by the factory of this class.
  */
  typedef ECELL_MODL *(*make_ecell_func) (GCELL *gcell);

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
  static ECELL_MODL *make_ecell (GCELL *gcell, string implementation);

 protected:
  
  /**
     Call this function to register a make function for the ecell
     based on the gcell type specified, and of given implementation.
  */
  static bool register_make_func (make_ecell_func make, string gcell_type, string implementation);

 public: // object functions ////////////////////////////////////////

  /**
   */
  ECELL_MODL (GCELL *gcell) : ECELL(gcell) {
    _mat = 0; _dofmappers.clear ();
  }
  /**
     Assemble the stiffness matrix.
   */
  virtual bool assemble_stiffness_matrix (FIELD_VECTOR *u, EVS_OOOFS *evs, FIELD_VECTOR *geometry) = 0;
  /**
    Assemble the mass matrix
   */ 
  virtual bool assem_consistent_mass (FIELD_VECTOR *u,EVS_OOOFS *evs, FIELD_VECTOR *geometry, double mmult)=0;
  
  
  virtual bool assemble_geom_stiffness_matrix (FIELD_VECTOR *u,
                           EVS_OOOFS *evs, FIELD_VECTOR *geometry, double pressure)    = 0;

  virtual bool assemble_geom_stiffness_matrix (FIELD_VECTOR *u,
                           EVS_OOOFS *evs, FIELD_VECTOR *geometry) = 0;
  virtual  double strain_energy (FIELD_VECTOR *u, FIELD_VECTOR *geometry) = 0 ;
  virtual bool calc_stress (POINT &param_loc,FIELD_VECTOR *u, FIELD_VECTOR *geometry ,double stress[6], double strain[6],double lambda, double mu) = 0 ;
  virtual double mass (FIELD_VECTOR *geometry) = 0;
  /**
     Attach the evaluation cell to the protocol.
     This includes:
     (a) get all needed parameters from the database;
     (b) establish a binding of the degrees of freedom at the nodes
         interacting at the ecell to the EVS.
     Arguments: database, evs (Eigen Value Solver), field name (what name you wish
     the solver to know the field under), field.
   */
  void attach (class PROTO_MODL *proto, string ggpath);
  /**
     Return the vector of dofmappers.
   */
  vector <DOFMAPPER_VECTOR *> dofmappers () { return _dofmappers; }
  
  /**
   */
  virtual ~ECELL_MODL () {};

  /**
     Set the name of material this cell should use.
  */
  void set_mat (MAT_ELASTICITY *mat) { _mat = mat; }
  /**
     Get the material the ecell should use.
     */
  MAT_ELASTICITY * mat () const { return _mat; }

 private: // class objects  /////////////////////////////////////////

  static OFACTORY<ECELL_MODL, make_ecell_func, GCELL *> *_ecell_factory;

 private: // object data ////////////////////////////////////////////

  vector <DOFMAPPER_VECTOR *>    _dofmappers;
  MAT_ELASTICITY                *_mat;
  
};

#endif
