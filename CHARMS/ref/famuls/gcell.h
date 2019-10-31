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
#ifndef GCELL_H
# define GCELL_H

#include <list>
#include <vector>
#include <set>
#include "fen.h"
#include "ofactory.h"
#include "conn.h"
#include "elixir_interface.h"
#include "fixed_vector.h"

class BFUN;

class GCELL {

 public: // declarations ////////////////////////////////////////////

  typedef CONN_MANIFOLD_DIM GCELL_MANIFOLD_DIM;
  static const sizet     ID_NONE = 0;
  /**
     Prototype of a `make' functions.
     The specializations of GCELL must supply a 'make' function
     so that instances of them may be made the GCELL factory.  The GCELL
     class provides a class method, make_gcell(), that creates
     a GCELL specialization instance of the specified type and of given
     implementation.
  */
  typedef class GCELL *(*make_gcell_func) (std::vector<FEN *> fens);

 public: // class functions /////////////////////////////////////////

  /**
     Register a function to make a gcell of a particular type.
   */
  static bool register_make_func (make_gcell_func make, std::string gcell_type, std::string implementation);
  /**
     Make a gcell.  This function uses a registered `make' function to
     actually produce a gcell of the requisite type and implementation.
     The function may return NULL if the factory does not have an
     appropriate make function, or if the make function did not get
     the proper number of finite element nodes in the fens vector.
  */
  static GCELL *make_gcell (std::string type, std::string implementation, std::vector <FEN *> fens);
  
 public: // object functions ////////////////////////////////////////
  
  GCELL ();
  virtual ~GCELL () {}
  /**
     Type name.  This is the name that is passed on to the make_gcell() function.
  */
  virtual const std::string type_name () { return std::string ("gcell"); }
  
  /**
     Return the gcell group to which this gcell belongs.
     A gcell may belong to a single group only.
  */
  class GCELL_GROUP *gcell_group () { return _gcell_group; }
  /**
     Set the gcell group to which this gcell belongs.
     A gcell may belong to a single group only.
  */
  void set_gcell_group (class GCELL_GROUP *gcell_group) { _gcell_group = gcell_group; }
  /**
     What is the manifold dimension of the gcell?
     (Point=0, line=1, surface=2, or solid=3.) 
   */
  virtual GCELL_MANIFOLD_DIM manifold_dim () = 0;
  /**
     What is the characteristic volume of the cell in the reference coordinates?
  */
  double char_dim ();
  /**
     Return a pointer to the connectivity this gcell is based upon.
     Each gcell is based upon some connectivity (but connectivity
     does not necessarily have to be associated with a gcell: consider
     the edges of a triangular gcell, they are not gcells, but they
     are connectivities).  Each gcell specialization uses a particular
     connectivity to express its manifold dimension, number of connected
     nodes, and the refinement pattern.
  */
  virtual const CONN_BASE *conn () = 0;

  /**
     Evaluate the basis function set associated with the evaluation point.
     It is assumed that the gcell on which this methods gets invoked
     is at the top of the hierarchy.  The reason is that the basis
     function evaluation proceeds from the finest level (top) to the
     coarsest level (bottom).
     The standard finite element basis functions defined
     over the gcell specialization (for instance, the tri-linear functions
     over an 8-node hexahedron) need to be evaluated for each cell in the
     hierarchy.  The basis functions (and their derivatives) are
     computed in the natural coordinate system.  When the whole hierarchy
     had been traversed, the Jacobian transformation is evaluated
     and the derivatives are transformed into the given (Cartesian)
     coordinate system. 
   */
  virtual bool eval_bfun_set (class EVALPT *evalpt) = 0;

 public: // these methods implement the gcell refinement hierarchy ////

  /**
     At which level in the hierarchy is the cell?
     The coarsest level is zero (0).
  */
  sizet level ();
  /**
     Return the parent of this gcell.  May be null if the
     cell has no parent (meaning it is at the bottom).
  */
  virtual GCELL *parent () = 0;
  /**
     Get the i-th child of the cell.
  */
  virtual GCELL *child (sizet i) const = 0;
  /**
     Get the number of children.  If the gcell is not refined, zero should be
     returned; otherwise the number of child gcells is returned
     (depending on the type of the gcell).
  */
  virtual sizet nchildren () const { return 0; }
  /**
     Has the gcell been divided into children?  Note we're not asking if they
     *could* exist, rather whether they actually exist.
  */
  virtual bool divided () const { return false; }
  /**
     Physically divide the gcell by creating its children (and all the nodes
     these children connect on the higher level).  Needs to be defined
     by the implementation of the specialization.
  */
  virtual void divide (class REF_CTX *ref_ctx) = 0;
  /**
     Return the set of detail functions that are needed for the refinement
     of the function associated with the node on input.
     This function may be called *ONLY* if the gcell had been divided
     into children before.
  */
  // virtual std::set <FEN *> detail_set (FEN *fen_of_bfun_to_refine) = 0;
 virtual void detail_set (FEN *fen_of_bfun_to_refine,set <FEN*> *rf) = 0;
  /**
     Return the set of all refinement functions (i.e. detail functions plus
     the one private function) that are needed for the refinement
     of the function associated with the node on input.
     This function may be called *ONLY* if the gcell had been divided
     into children before.
  */
  // virtual std::set <FEN *> complete_refinement_set (FEN *fen_of_bfun_to_refine) = 0;
  virtual void  complete_refinement_set (FEN *fen_of_bfun_to_refine,set <FEN*> *rf) = 0;
  /**
     Map a parametric point xi, eta, theta (one, two, or all three of these
     coordinates may be used) into a child of this gcell.
  */
  virtual bool
    map_to_child (POINT & param_loc,
                  GCELL **child,
                  POINT *child_param_loc) = 0;
  /**
    Given a parametric location, return the parent and the parametric coordinates
     to which this location maps in the parent of the this cell.
     If the cell has no parent, false is returned and the output parameters are undefined.
     If successful, *parent is the parent of the cell, and
     true is returned.
  */ 
  virtual bool
    map_to_parent (POINT & param_loc,
                   GCELL **parent, POINT *parent_param_loc) = 0;
  /**
     Map a finite element node to the parametric coordinate in the
     gcell given as argument.  If the finite element node happens
     not to be referenced by the connectivity of the cells, an
     exception EXCEPTION_BAD_VALUE is thrown.
  */
  virtual void map_fen (FEN *fen, POINT *param_loc) = 0;
  /**
     Is this the topmost gcell that is used by any basis function in the field?
  */
  bool is_topmost_gcell_in_field (class FIELD_BASE *field);
  /**
     Find the topmost gcell that is used by any basis function in the
     field given as argument.  If no such gcell exists, null (0) is returned.
     If non-null is returned, *txi, *teta, and *ttheta return the parametric
     coordinates of the point xi, eta, theta mapped into the topmost gcell.
     Note: This is a generic procedure that relies on the implementations
     of map_to_child() and map_to_parent() in the specializations of the gcell.
  */
  GCELL * map_to_topmost_gcell_in_field (class FIELD_BASE *field,
                                         POINT &param_loc,
                                         POINT *mapped_to_param_loc);
  /**
     Find the root gcell of the refinement hierarchy whose member `this' is.
     (It's the Adam or Eve ;)
  */
  GCELL *root() {
    GCELL *r = this;
    while (r->parent()) {
      r = r->parent();
    }
    return r;
  }


 public: // debugging methods ///////////////////////////////////////

  /**
     Print out the gcell definition.  No newline.
   */
  virtual void debug_display ();
#if defined (GIFACE) && GIFACE
  /**
     Create an Elixir graphics object.
   */
  virtual GraphicObj *make_elixir_obj () { return 0; }
  virtual GraphicObj *make_elixir_obj (vector <FIXED_VECTOR<3> > &p) { return 0; }
#endif

 void set_id(sizet id) {
    CHECK (_id == ID_NONE,EXCEPTION_ILLEGAL_USE,;);
    _id = id;
  }

 sizet id() { return _id;}
 private: // class data /////////////////////////////////////////////

  static OFACTORY<GCELL, make_gcell_func, std::vector <FEN *> > *_gcell_factory;

 private: // object data ////////////////////////////////////////////

  class GCELL_GROUP     *_gcell_group;
  sizet                  _id;
  
};

#endif
