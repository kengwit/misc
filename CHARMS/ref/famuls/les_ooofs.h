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
#ifndef LES_OOOFS_H
# define LES_OOOFS_H

#include <map>
extern "C" {
#include "fles.h"
}
#include "dofmapper_scalar.h"
#include "dofmapper_vector.h"
#include "field.h"
#include "mgr.h"
#include "ooofs/ooofs.h"


class LES_OOOFS {

 public: // object functions ////////////////////////////////////////

  /**
     Create a linear equation solver (LES).  Initialize from the database.
   */
  LES_OOOFS (string name, MGR *mgr, string options = "");

  ~LES_OOOFS () {
    if (_K) delete _K;
    for (dofmappers_map::const_iterator i = _dofmappers_map.begin ();
         i != _dofmappers_map.end(); i++) {
      dofmapper_map *m = i->second;
      for (dofmapper_map::const_iterator j = m->begin (); j != m->end(); j++) {
        DOFMAPPER_BASE *dofmapper = j->second;
        delete dofmapper;
      }
      delete m;
    }
    if (_col_heights) delete [] _col_heights;
  }
  /**
     Return the total number of equations that this linear equation solver
     presently handles.
   */
  int tot_of_equations ();
  /**
     Start assembly by zeroing out the lhs matrix.
  */
  void zero_lhs ();
  /**
     Start assembly by zeroing out the rhs vector.
  */
  void zero_rhs ();
    /**
     Add a submatrix to the global matrix.  The submatrix couples
     the degrees of freedom given as the first two arguments.
  */
  template <int NUM_COMPONENTS1, int NUM_COMPONENTS2>
    void lhs_add (DOFMAPPER<NUM_COMPONENTS1> *dofmap1, DOFMAPPER<NUM_COMPONENTS2> *dofmap2,
                       double values[NUM_COMPONENTS1][NUM_COMPONENTS2])
    {
      CHECK (_K != NULL, EXCEPTION_NULL_PTR,;);
      
      if (_total_of_equations <= 0) tot_of_equations ();
      int rows[NUM_COMPONENTS1], cols[NUM_COMPONENTS2];
      dofmap1->eqnums (rows);
      dofmap2->eqnums (cols);
      
      OOOFS::IntArray ia(2);
      OOOFS::FloatMatrix fm(2,2);  fm.zero();
      
      for (int r = 0; r < NUM_COMPONENTS1; r++) {
        if (rows[r] != DOFMAPPER_BASE::invalid_eqnum) {
          for (int c = 0; c < NUM_COMPONENTS2; c++) {
            if (cols[c] != DOFMAPPER_BASE::invalid_eqnum) {
              ia(0)=rows[r]; ia(1)=cols[c];
              fm(0,1)=values[r][c];
              _K->assemble (ia, fm);
            }
          }
        }
      }
    }
      
  /**
     Assemble a subvector to the global rhs vector.
  */
  template <int NUM_COMPONENTS>
    void rhs_add (DOFMAPPER<NUM_COMPONENTS> *dofmap, double values[NUM_COMPONENTS])
    {
      if (_total_of_equations <= 0) tot_of_equations ();
      int rows[NUM_COMPONENTS];
      dofmap->eqnums (rows);
      // now check that all equation numbers are valid; if they are not
      // delete the row (column) and assemble those zeros (arbitrarily)
      // to equation zero.
      OOOFS::IntArray ia(1);
      OOOFS::FloatArray fa(1);  fa.zero();
      for (int r = 0; r < NUM_COMPONENTS; r++) {
        if (rows[r] != DOFMAPPER_BASE::invalid_eqnum) {
          ia.at(1)=rows[r];
          fa.at(1)=values[r];
          _F.assemble(fa, ia);
        }
      }
    }
  /**
     Finish the assembly of the rhs vector.
  */
  void finish_rhs ();
  /**
     Finish the assembly of the lhs matrix.
  */
  void finish_lhs ();
  /**
     Solve the global system.
  */
  bool solve ();
  /**
   */
  void dump_lhs (string file) { }
  /**
   */
  void dump_rhs (string file) {  }
  /**
   */
  void dump_sol (string file) {  } 
  
  /**
     Distribute the solution of the global system to the
     field on input.  The field_name parameter refers to the
     name under which the field is known to the solver.
     It doesn't have to have anything to do with the name
     under which this field is known to the rest of the program.
     Field_name has to be the same field_name that gets passed to
     LES::dofmapper().
  */
  template <int NUM_COMPONENTS>
    bool solution (string field_name, FIELD<NUM_COMPONENTS> *field)
    {
      // first find the dofmapper map
      dofmapper_map *m;
      {
        dofmappers_map::const_iterator i = _dofmappers_map.find (field_name);
        CHECK (i != _dofmappers_map.end (), EXCEPTION_BAD_VALUE,;);
        m = i->second;
      }
      // now distribute 
      dofmapper_map::const_iterator i = m->begin ();
      while (i != m->end ()) {
        const BFUN_DOFPARAM_PAIR_ID dofparam_id = i->first;
        const DOFMAPPER_BASE *dofmapper = i->second;
        const int neqns = dofmapper->neqns ();
        int *eqnums = new int [neqns];
        dofmapper->eqnums (eqnums);
        for (int j = 0; j < neqns; j++) {
          if (eqnums[j] != DOFMAPPER_BASE::invalid_eqnum) {
            double value = _F.at(eqnums[j]);
            field->field_pair (dofparam_id)->set_dofparam_comp (j, value);
          }
        }
        delete [] eqnums;
        i++;
      }
      return true; // ok
    }

  /**
     Generate a dofmapper object for the named field and given dofparam_id.
     This creates a mapping object that an ecell may store without
     needing to know the global numbering of the dofs.
     The LES knows the names of the fields it is constructed to deal with,
     but it has no idea what those fields are.  Field_name has to be
     the same field_name that gets passed to LES::solution().
     

     First, the LES checks if it already has a dofmapper for the dofparam_id
     given on input.  If it does, the existing dofmapper is returned;
     otherwise, a new one is generated by the field, and returned to the caller.

     Note: the receiver of the dofmapper must not delete the dofmapper.
     These are shared objects (many clients may be given the same dofmapper),
     and it is therefore the task of the LES to delete the dofmappers
     when the LES itself is destroyed.
  */
  template <int NUM_COMPONENTS>
    DOFMAPPER<NUM_COMPONENTS> *dofmapper (string field_name,
                                               FIELD<NUM_COMPONENTS> *field,
                                               BFUN_DOFPARAM_PAIR_ID dofparam_id) {
    // first find the dofmapper map
    dofmapper_map *m;
    {
      dofmappers_map::const_iterator i = _dofmappers_map.find (field_name);
      if (i == _dofmappers_map.end ()) { // field not known to the LES yet
        m = new dofmapper_map;
        _dofmappers_map.insert (dofmappers_map::value_type (field_name, m));
      } else
        m = i->second;
    }
    // now produce the dofmapper
    dofmapper_map::iterator i = m->find (dofparam_id);
    DOFMAPPER<NUM_COMPONENTS> *dofmapper = 0;
    if (i == m->end ()) {   // not in the map yet
      dofmapper = field->dofmapper (dofparam_id);
      m->insert (dofmapper_map::value_type (dofparam_id, dofmapper));  // add it
    } else {
      dofmapper = dynamic_cast <DOFMAPPER<NUM_COMPONENTS> *> (i->second);
    }
    if (_making_conn_graph) add_to_batch (dofmapper);
    return dofmapper;
  }
  
  /**
    Return the name of this object.
  */
  const string name () { return _name; }
  /**
     Start gathering the new batch of connections among basis functions.
  */
  void conn_graph_new_batch ();
  /**
     Update the connection graph.
  */
  void update_conn_graph ();
  
  
 private: // object data ////////////////////////////////////////////

  typedef map <BFUN_DOFPARAM_PAIR_ID, DOFMAPPER_BASE *> dofmapper_map;
  typedef map <string, dofmapper_map *> dofmappers_map;
  typedef map <DOFMAPPER_BASE *, set <DOFMAPPER_BASE *> > conn_graph_map;
  
  string                      _name;
  MGR                        *_mgr;
  string                      _options;
  dofmappers_map              _dofmappers_map;
  int                         _total_of_equations;
  OOOFS::SparseMtrx          *_K;
  OOOFS::FloatArray           _F;
  conn_graph_map              _conn_graph_map;
  bool                        _making_conn_graph;
  set <DOFMAPPER_BASE *>      _batch;
  int                        *_col_heights;

 private: // object functions ///////////////////////////////////////

  void setup_solver_context (DB *db, string options);
  void naive_equation_numbering ();
  void make_col_heights ();
  void add_to_batch (DOFMAPPER_BASE *dm);

  class DOFMAPPER_BASE_COMPARATOR {
  public:
    bool operator() (DOFMAPPER_BASE *p0, DOFMAPPER_BASE *p1) const {
      return (p0->fen()->id() < p1->fen()->id());
    }
  };

};

#endif
