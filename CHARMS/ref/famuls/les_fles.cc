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
#include "famuls.h"
#include "args.h"
#include "les_fles.h"

static int _tot_of_equations (void *client)
{
  LES_FLES *a = (LES_FLES *)client;
  int n = a->tot_of_equations (); 
  return n;
}

static int _error_notify (fles_error_enum error_num,
                          int line,
                          char *file,
                          char *detail,
                          void *client_data)
{
  char buf[1024];
  sprintf (buf, "%s (l.%d), %s", file, line, detail);
  CHECK (0, EXCEPTION_IN_LOCAL_LIBS, (buf));
  return 0;
}

void
LES_FLES::setup_solver_context (DB *db, string options)
{
  string fles_solver = "gep"; // default: basic Gauss elimination solver
  double nz_per_row  = 0.1; // hierarchical basis leads to denser matrices
  bool use_nnzs = false;

  if (db) {
    string path = "les/" + this->name () + "/fles_solver";
    fles_solver = db->DB_GET_STRING (path);
  }
  
  if (!options.empty ()) {
    tokens_parser parser = tokens_new_parser ();
    CHECK (parser, EXCEPTION_NULL_PTR,;);
    const char *opts = options.c_str();
    if (tokens_parse_line (parser, (char *) opts)) {
      int pos = 0;
      if ((pos = tokens_have_keyword (parser, 1, "NZ_PER_ROW")) > 0) {
        if (tokens_total_of_tokens (parser) >= pos+1) {
          nz_per_row = tokens_token_as_double (parser, pos+1);
        }
      }
      if ((pos = tokens_have_keyword (parser, 1, "USE_NNZS")) > 0) {
        use_nnzs = true;
      }
      if ((pos = tokens_have_keyword (parser, 1, "FLES_SOLVER")) > 0) {
        if (tokens_total_of_tokens (parser) >= pos+1) {
          const char *solvname;
          solvname = tokens_token_as_string (parser, pos+1);
          if        (tokens_token_and_keyword_equiv (solvname, "SMESCHACH")) {
            fles_solver = "smeschach";
          } else if (tokens_token_and_keyword_equiv (solvname, "DMESCHACH")) {
            fles_solver = "dmeschach";
#         if defined (SOLVER_PETSC) && SOLVER_PETSC
          } else if (tokens_token_and_keyword_equiv (solvname, "SPETSC")) {
            fles_solver = "spetsc";
#         endif
          } else if (tokens_token_and_keyword_equiv (solvname, "GEP")) {
            fles_solver = "gep";
          }
        }
      }
    }
    tokens_delete_parser (parser);
  } // have options

  if (fles_solver == "smeschach") {
    int method = FLES_SMESCHACH_SPCG;
    _slesctx = fles_new_sparse_meschach_solver (this, _tot_of_equations, _error_notify);
    _slesctx->set_option (_slesctx, FLES_SMESCHACH_METHOD, (void *)method);
#if defined(SOLVER_PETSC) && SOLVER_PETSC
  } else if (fles_solver == "spetsc") {
    _slesctx = fles_new_spetsc_solver (this, _tot_of_equations, _error_notify);
    _slesctx->set_option (_slesctx, PETSC_OPTS_NZ_PER_ROW, (void *)&nz_per_row);
#endif
  } else if (fles_solver == "dmeschach") {
    _slesctx = fles_new_dense_meschach_solver (this, _tot_of_equations, _error_notify);
  } else if (fles_solver == "gep") {
    _slesctx = fles_new_gep_solver (this, _tot_of_equations, _error_notify);
    _slesctx->set_option (_slesctx, GEP_OPTS_RCOND, (void *)0);
  } else {
    NOT_IMPLEMENTED ("unknown FLES solver");
  }

  CHECK (_slesctx, EXCEPTION_NULL_PTR,;);
}

LES_FLES::LES_FLES (string name, MGR *mgr, string options) 
{
  _name = name;
  _mgr  = mgr;
  _options = options;
  _slesctx = 0;
  _dofmappers_map.clear ();
  _total_of_equations = 0;
  _conn_graph_map.clear ();
  _making_conn_graph = false;
  _nnzs = 0;
}


void LES_FLES::zero_lhs () {
  if (_slesctx == 0) setup_solver_context (_mgr->db(), _options);
  _slesctx->zero_lhs (_slesctx);
} 

void LES_FLES::zero_rhs () {
  if (_slesctx == 0) setup_solver_context (_mgr->db(), _options);
  _slesctx->zero_rhs (_slesctx);
} 

void LES_FLES::finish_lhs ()
{
  CHECK (_slesctx->finish_lhs != NULL, EXCEPTION_NULL_PTR,;);
  _slesctx->finish_lhs (_slesctx, FLES_FINAL_FINISH);
  // _slesctx->dump_lhs (_slesctx, "lhs");
}

void LES_FLES::finish_rhs ()
{
  CHECK (_slesctx->finish_rhs != NULL, EXCEPTION_NULL_PTR,;);
  _slesctx->finish_rhs (_slesctx);
  // _slesctx->dump_rhs (_slesctx, "rhs");
}

void LES_FLES::naive_equation_numbering ()
{
  _total_of_equations = 0;
  dofmappers_map::const_iterator fi = _dofmappers_map.begin ();
  while (fi != _dofmappers_map.end ()) {
    dofmapper_map *m = fi->second;
    dofmapper_map::const_iterator i = m->begin ();
    while (i != m->end ()) {
      DOFMAPPER_BASE *dofmapper = i->second;
      const int neqns = dofmapper->neqns ();
      int *en = new int [neqns];
      bool *c = new bool [neqns];
      // first get the current equation constraints
      dofmapper->constrained (c);
      // now number the unconstrained equations
      for (int j = 0; j < neqns; j++) {
        if (!c[j]) { en[j] = _total_of_equations; _total_of_equations++; }
        else       { en[j] = DOFMAPPER_BASE::invalid_eqnum; }
      }
      // and set the equation numbers in the dofmapper
      dofmapper->set_eqnums (en);
      delete [] en;
      delete [] c;
      i++;
    }
    fi++;
  }
  if (! _conn_graph_map.empty()) make_nnzs ();
}

int LES_FLES::tot_of_equations ()
{
  if (_total_of_equations <= 0) {
    naive_equation_numbering ();
  }
  CHECK (_total_of_equations > 0, EXCEPTION_BAD_VALUE,;);
  return _total_of_equations;
}

bool LES_FLES::solve ()
{
  //string file_name = (char*) _ncycle ;//RAW :ashi
 // write_lhs(_slesctx,"k.jpp" ); //RAW :ashi
 // _ncycle ++; //RAW :ashi
  CHECK (_slesctx->solve, EXCEPTION_NULL_PTR,;);
  if (!_slesctx->solve (_slesctx)) {
    // RAW: recovery?
    return false;
  }
  return true;
}

void LES_FLES::conn_graph_new_batch ()
{
  _making_conn_graph = true;
  _batch.clear();
}

void LES_FLES::update_conn_graph ()
{
  for (set <DOFMAPPER_BASE *>::iterator i = _batch.begin();
       i != _batch.end(); i++) {
    conn_graph_map::iterator cgmi = _conn_graph_map.find (*i);
    set <DOFMAPPER_BASE *> s;
    if (cgmi == _conn_graph_map.end()) {
      set <DOFMAPPER_BASE *> emptyset; emptyset.clear();
      _conn_graph_map.insert (conn_graph_map::value_type (*i, emptyset));
      cgmi = _conn_graph_map.find (*i);
    }
    for (set <DOFMAPPER_BASE *>::iterator j = _batch.begin();
         j != _batch.end(); j++) {
      (*cgmi).second.insert (*j);
    }
  }
  _making_conn_graph = false;
}

void LES_FLES::add_to_batch (DOFMAPPER_BASE *dm)
{
  _batch.insert (dm);
}

void LES_FLES::make_nnzs ()
{
  // Allocate and initialize the array
  _nnzs = new int [_total_of_equations];
  for (int k = 0; k < _total_of_equations; k++) _nnzs[k] = 0;

  // Now loop over all dofmappers
  dofmappers_map::const_iterator fi = _dofmappers_map.begin ();
  while (fi != _dofmappers_map.end ()) {
    dofmapper_map *m = fi->second;
    dofmapper_map::const_iterator i = m->begin ();
    while (i != m->end ()) {
      DOFMAPPER_BASE *dofmapper = i->second;
      const int neqns = dofmapper->neqns();
      int *eqnums = new int [neqns];
      dofmapper->eqnums (eqnums);
      conn_graph_map::iterator cgmi = _conn_graph_map.find(dofmapper);
      CHECK (cgmi != _conn_graph_map.end(), EXCEPTION_NULL_PTR,;);
      set <DOFMAPPER_BASE *> &s = (*cgmi).second;
      for (set <DOFMAPPER_BASE *>::iterator si = s.begin(); si != s.end(); si++) {
        DOFMAPPER_BASE *dm = *si;
        const int neqns2 = dofmapper->neqns();
        int * eqnums2 = new int [neqns2];
        dm->eqnums (eqnums2);
        for (int k = 0; k < neqns; k++) {
          if (eqnums[k] != DOFMAPPER_BASE::invalid_eqnum) {
            for (int j = 0; j < neqns2; j++) {
              if (eqnums2[j] != DOFMAPPER_BASE::invalid_eqnum) {
                _nnzs[eqnums[k]] += 1;
              }
            }
          }
        }      
        delete [] eqnums2;
      }
      delete [] eqnums;
      i++;
    }
    fi++;
  }

//   cerr << "Conn" << endl;
//   cerr << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
//   for (int k = 0; k < _total_of_equations; k++) {
//     cerr << k << ":" <<_nnzs[k] << endl;
//   }
//   cerr << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
//   cerr << "use_nnzs set (" << _nnzs << ") (" << __FILE__ << ")" << endl;

#if defined (SOLVER_PETSC) && SOLVER_PETSC
  _slesctx->set_option (_slesctx, PETSC_OPTS_NNZ_PER_ROW, (void *)_nnzs);
#endif
}
