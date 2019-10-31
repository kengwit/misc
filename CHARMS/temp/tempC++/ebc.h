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
#ifndef EBC_H
# define EBC_H

#include "field.h"
#include "selcrit_fen.h"
#include "selcrit_box.h"
#include "selcrit_line.h"
#include "selcrit_quad.h"
#include "selcrit_ntri.h"
#include "selcrit_tri.h"
#include <string>
#include <list>
#include <map>
extern "C" {
#include "expevalP.h"
#include "tokensP.h"
}
extern "C" {
#include "lobs3d.h"
}
#include "func.h"

template <int NUM_COMPONENTS>
class EBC {

public:

  EBC (GMESH *gmesh) {
    _gmesh = gmesh;
    _selcrits.clear();
    _lobs = lobs_3d_new_instance (0, 5);
    CHECK (_lobs, EXCEPTION_NULL_PTR,;);
  }

  ~EBC () {
    for (typename list <FUNC_SELCRIT_PAIR *>::iterator i = _selcrits.begin();
      i != _selcrits.end(); i++) {
        delete (*i);
      }
      if (_lobs) lobs_3d_free_instance (_lobs);
  }
  /**
  Apply EBC to the field given as argument.
  The EBC object maintains a search structure (LOBS: Library for
  Overlapping Object Searching) to speed up application of the
  boundary conditions.  Each selection criterion is associated
  with a bounding box  which is stored in the search structure.
  Given the location of a node, all potentially applicable criteria
  may be found quickly using the search structure.
  */
  void apply (FIELD<NUM_COMPONENTS> *field);
  /**
  Apply time-dependent EBC to the field given as argument.
  The EBC object maintains a search structure (LOBS: Library for
  Overlapping Object Searching) to speed up application of the
  boundary conditions.  Each selection criterion is associated
  with a bounding box  which is stored in the search structure.
  Given the location of a node, all potentially applicable criteria
  may be found quickly using the search structure.
  */
  void apply (FIELD<NUM_COMPONENTS> *field, double time);
  /**
  Read the EBC definition from a file.
  */
  bool read (string file);

private: // object data ////////////////////////////////////////////

  class FUNC_SELCRIT_PAIR {
  public:
    /**
    The function/selection-criterion pair *owns* the
    pointers!  It will delete both in the destructor.
    */
    FUNC_SELCRIT_PAIR (FUNC<NUM_COMPONENTS> *func, SELCRIT_BASE *selcrit) : _func(func), _selcrit(selcrit) {}
    ~FUNC_SELCRIT_PAIR () {
      delete this->_func;
      delete this->_selcrit;
    }
    FUNC<NUM_COMPONENTS> *_func;
    SELCRIT_BASE         *_selcrit;
  };

  GMESH                     *_gmesh;
  list <FUNC_SELCRIT_PAIR *> _selcrits;
  lobs_3d_front_end          _lobs;

private: // object functions ///////////////////////////////////////

  void add (FUNC_SELCRIT_PAIR *fp) {
    _selcrits.push_back(fp);
    G3D_box_t b = fp->_selcrit->box (_gmesh);
    lobs_3d_add_object (_lobs, fp, b.lower.x, b.lower.y, b.lower.z,
      b.upper.x, b.upper.y, b.upper.z);
  }

};


template <int NUM_COMPONENTS>
bool
EBC<NUM_COMPONENTS>::read (string file)
{
  tokens_parser parser = tokens_new_parser ();
  CHECK (parser, EXCEPTION_NULL_PTR,;);

  CHECK (tokens_open_file (parser, (char *) file.c_str ()), EXCEPTION_FILE_ERROR, (file));

  // Initialize the expressions
  string expr[NUM_COMPONENTS];
  for (sizet j = 0; j < NUM_COMPONENTS; j++) {
    expr[j] = "";
  }

  while (tokens_next_line (parser)) {
    CHECK (tokens_total_of_tokens (parser) >= 1, EXCEPTION_BAD_VALUE,;);
    if (tokens_token_and_keyword_equiv (tokens_token_as_string (parser, 1), "function")) {
      CHECK (tokens_total_of_tokens (parser) == NUM_COMPONENTS+1, EXCEPTION_BAD_VALUE,;);
      for (sizet j = 0; j < NUM_COMPONENTS; j++) {
        expr[j] = string (tokens_token_as_string (parser, 2+j));
      }
    } else if (tokens_token_and_keyword_equiv (tokens_token_as_string (parser, 1), "fen")) {
      CHECK (tokens_total_of_tokens (parser) == 2, EXCEPTION_BAD_VALUE,;);
      int id = tokens_token_as_int (parser, 2);
      FUNC_SELCRIT_PAIR *fp
        = new FUNC_SELCRIT_PAIR (new FUNC<NUM_COMPONENTS>(expr), new SELCRIT_FEN(id));
      add (fp);
    } else if (tokens_token_and_keyword_equiv (tokens_token_as_string (parser, 1), "box")) {
      CHECK (tokens_total_of_tokens (parser) == 9, EXCEPTION_BAD_VALUE,;);
      POINT p;
      p(0) = tokens_token_as_double (parser, 2);
      p(1) = tokens_token_as_double (parser, 3);
      p(2) = tokens_token_as_double (parser, 4);
      POINT q;
      q(0) = tokens_token_as_double (parser, 5);
      q(1) = tokens_token_as_double (parser, 6);
      q(2) = tokens_token_as_double (parser, 7);
      double inflate = tokens_token_as_double (parser, 9);
      FUNC_SELCRIT_PAIR *fp
        = new FUNC_SELCRIT_PAIR (new FUNC<NUM_COMPONENTS>(expr), new SELCRIT_BOX (p, q, inflate));
      add (fp);
    } else if (tokens_token_and_keyword_equiv (tokens_token_as_string (parser, 1), "quad")) {
      // quad corner 1 ... corner 4        dist  d function f
      //  1   2 3 4 5 6 7 8 9 10 11 12 13   14  15  16      17
      CHECK (tokens_total_of_tokens (parser) == 15, EXCEPTION_BAD_VALUE,;);
      POINT p, q, r, s;
      p(0) = tokens_token_as_double (parser, 2);
      p(1) = tokens_token_as_double (parser, 3);
      p(2) = tokens_token_as_double (parser, 4);
      q(0) = tokens_token_as_double (parser, 5);
      q(1) = tokens_token_as_double (parser, 6);
      q(2) = tokens_token_as_double (parser, 7);
      r(0) = tokens_token_as_double (parser, 8);
      r(1) = tokens_token_as_double (parser, 9);
      r(2) = tokens_token_as_double (parser, 10);
      s(0) = tokens_token_as_double (parser, 11);
      s(1) = tokens_token_as_double (parser, 12);
      s(2) = tokens_token_as_double (parser, 13);
      double dist = tokens_token_as_double (parser, 15);
      FUNC_SELCRIT_PAIR *fp
        = new FUNC_SELCRIT_PAIR (new FUNC<NUM_COMPONENTS>(expr),
        new SELCRIT_QUAD (p, q, r, s, dist));
      add (fp);
    } else if (tokens_token_and_keyword_equiv (tokens_token_as_string (parser, 1), "tri")) {
      // tri corner 1 ... corner 3         dist  d 
      //  1   2 3 4 5 6 7 8 9 10            11   12 
      CHECK (tokens_total_of_tokens (parser) == 12, EXCEPTION_BAD_VALUE,;);
      POINT p, q, r;
      p(0) = tokens_token_as_double (parser, 2);
      p(1) = tokens_token_as_double (parser, 3);
      p(2) = tokens_token_as_double (parser, 4);
      q(0) = tokens_token_as_double (parser, 5);
      q(1) = tokens_token_as_double (parser, 6);
      q(2) = tokens_token_as_double (parser, 7);
      r(0) = tokens_token_as_double (parser, 8);
      r(1) = tokens_token_as_double (parser, 9);
      r(2) = tokens_token_as_double (parser, 10);
      double dist = tokens_token_as_double (parser, 12);
      FUNC_SELCRIT_PAIR *fp
        = new FUNC_SELCRIT_PAIR (new FUNC<NUM_COMPONENTS>(expr), new SELCRIT_TRI (p, q, r, dist));
      add (fp);
    } else if (tokens_token_and_keyword_equiv (tokens_token_as_string (parser, 1), "ntri")) {
      // ntri n0 n1 n2  dist  d
      //  1   2 3 4 5 6 
      CHECK (tokens_total_of_tokens (parser) == 6, EXCEPTION_BAD_VALUE,;);
      int n0, n1, n2;
      n0 = tokens_token_as_int (parser, 2);
      n1 = tokens_token_as_int (parser, 3);
      n2 = tokens_token_as_int (parser, 4);
      double dist = tokens_token_as_double (parser, 6);
      FUNC_SELCRIT_PAIR *fp
        = new FUNC_SELCRIT_PAIR (new FUNC<NUM_COMPONENTS>(expr),
        new SELCRIT_NTRI (n0, n1, n2, dist));
      add (fp);
    } else if (tokens_token_and_keyword_equiv (tokens_token_as_string (parser, 1), "line")) {
      CHECK (tokens_total_of_tokens (parser) == 9, EXCEPTION_BAD_VALUE,;);
      POINT p, q;
      p(0) = tokens_token_as_double (parser, 2);
      p(1) = tokens_token_as_double (parser, 3);
      p(2) = tokens_token_as_double (parser, 4);
      q(0) = tokens_token_as_double (parser, 5);
      q(1) = tokens_token_as_double (parser, 6);
      q(2) = tokens_token_as_double (parser, 7);
      double dist = tokens_token_as_double (parser, 9);
      FUNC_SELCRIT_PAIR *fp
        = new FUNC_SELCRIT_PAIR (new FUNC<NUM_COMPONENTS>(expr), new SELCRIT_LINE (p, q, dist));
      add (fp);
    }
  }
  tokens_close_curr_file (parser);
  tokens_delete_parser (parser);

  return true;
}

#include "evalpt.h"

template <int NUM_COMPONENTS>
void
EBC<NUM_COMPONENTS>::apply (FIELD<NUM_COMPONENTS> *field, double time)  {
  const double eps = 1000000 * 2e-16; // RAW entirely arbitrary 
  sizet maxlevel = field->max_level (field) ;
  const sizet npairs = field->npairs();

  for (sizet level = 0; level < maxlevel+1; level++) {

    for (sizet j=0; j<npairs; j++) {
      FIELD_PAIR<NUM_COMPONENTS> *fp = field->jth_field_pair (j);
      BFUN *bfun = fp->bfun();
      if (bfun->level() == level) {
        // set the constraints and real values on the boundary nodes
        FEN *fen = bfun->fen();
        POINT ref_loc = fen->ref_loc();
        lobs_3d_search (_lobs,
          ref_loc(0)-eps, ref_loc(1)-eps, ref_loc(2)-eps,
          ref_loc(0)+eps, ref_loc(1)+eps, ref_loc(2)+eps);
        FUNC_SELCRIT_PAIR *funpair
          = static_cast <FUNC_SELCRIT_PAIR *> (lobs_3d_first_overlapping_object (_lobs));
        while (funpair) { // Now loop over all  selection criteria, try to match
          if (funpair->_selcrit->match (fen)) {
            BFUN_DOFPARAM_PAIR_ID dpid = field->dofparam_id (fen);
            CHECK (dpid != INVALID_BFUN_DOFPARAM_PAIR_ID, EXCEPTION_BAD_VALUE,;);
            FIELD_PAIR<NUM_COMPONENTS> *fp = field->field_pair (dpid);
            POINT ref_loc = fen->ref_loc ();
            POINT4 at4; at4(0) = ref_loc(0); at4(1) = ref_loc(1); at4(2) = ref_loc(2); at4(3) = time;
            FIXED_VECTOR<NUM_COMPONENTS> v = (*funpair->_func) (at4);
            for (sizet j=0; j < NUM_COMPONENTS; j++) { // Set constraint tags
              if (funpair->_func->given(j)) {
                fp->set_constrained_comp (j, true, v(j));
              }
            }
            FIXED_VECTOR<NUM_COMPONENTS> srcv = fp->dofparam ();
            //Now get the nodal parameters for the boundary nodes.
            GCELL *gcell = bfun->gcell(0); // any gcell would do
            // Evaluate basis function values at level lower than the current level
            EVALPT evalpt (field, gcell, fen);
            evalpt.eval ();
            // Loop over all basis functions in the field at evalpt
            double N_j = 1;  // Adjust for PU 
            sizet nbfuns = evalpt.nbfuns();
            for (sizet j=0; j<nbfuns; j++) {
              BFUN *dbfun = evalpt.bfun(j);
              if (dbfun != bfun) {
                double N = evalpt.N(j);
                BFUN_DOFPARAM_PAIR_ID dpid = field->dofparam_id (dbfun->fen());
                CHECK (dpid != INVALID_BFUN_DOFPARAM_PAIR_ID, EXCEPTION_NULL_PTR,;);
                FIXED_VECTOR<NUM_COMPONENTS> dp = field->field_pair(dpid)->dofparam();
                srcv.add(-N, dp);
              } else {
                N_j = evalpt.N(j);
              }// if (dbfun
            }// for (sizet j=0; j<nbfuns; j++)
            srcv.scale(1/N_j); // Adjust for PU 
            fp->set_dofparam_constrained_only (srcv); // Set the constrained values
          }// if (funpair->_selcrit->match (fen))
          funpair = static_cast <FUNC_SELCRIT_PAIR *> (lobs_3d_next_overlapping_object (_lobs));
        }// while (funpair)

      }//if (bfun->level() == level)
    }// for (sizet j=0; j<npairs; j++)

  }// for (sizet level = 0; level < maxlevel+1; level++)
}

template <int NUM_COMPONENTS>
void
EBC<NUM_COMPONENTS>::apply (FIELD<NUM_COMPONENTS> *field)  {
  apply(field,0);
}

#endif
