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
#include "evs_ooofs.h"

static int _tot_of_equations (void *client)
{
  EVS_OOOFS *a = (EVS_OOOFS *)client;
  int n = a->tot_of_equations (); 
  return n;
}

static int _error_notify (int error_num,
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

double
EVS_OOOFS::setup_solver_context (DB *db, string options)
{
  return 0;  
}

EVS_OOOFS::EVS_OOOFS (string name, MGR *mgr, string options) 
{
  _name = name;
  _mgr   = mgr;
  _options = options;
  _K = 0;
  _M = 0;
  _add_to_K = true;
  _numberOfRequiredEigenValues = 0;
  _rtolv = 1e-3;
  _dofmappers_map.clear ();
  _total_of_equations = 0;
  _conn_graph_map.clear ();
  _making_conn_graph = false;
  _col_heights = 0;
  OOOFS::error_notify = _error_notify;
}

void EVS_OOOFS::zero_k () {
  if (_K == 0) {
    _K = new OOOFS::Skyline ();
    CHECK (_K, EXCEPTION_NULL_PTR,;);
    CHECK(tot_of_equations() > 0, EXCEPTION_BAD_VALUE,;);
    SMART_HANDLE<LOGGER_STREAM> lsb = _mgr->logger("building K", true);
    _K->buildInternalStructure (*lsb, _total_of_equations, _col_heights);
    *lsb << "done" << endl;
  }
  _K->zero();
} 

void EVS_OOOFS::zero_m () {
  if (_M == 0) {
    _M = new OOOFS::Skyline ();
    CHECK (_M, EXCEPTION_NULL_PTR,;);
    CHECK(tot_of_equations() > 0, EXCEPTION_BAD_VALUE,;);
    SMART_HANDLE<LOGGER_STREAM> lsb = _mgr->logger("building M", true);
    _M->buildInternalStructure (*lsb, _total_of_equations, _col_heights);
    *lsb << "done" << endl;
  }
  _M->zero();
} 

void EVS_OOOFS::finish_k ()
{
}

void EVS_OOOFS::finish_m ()
{
}

#if 1
void EVS_OOOFS::naive_equation_numbering ()
{
  _total_of_equations = 0;
  sizet n = 0;
  dofmappers_map::const_iterator fi = _dofmappers_map.begin ();
  while (fi != _dofmappers_map.end ()) {
    dofmapper_map *m = fi->second;
    dofmapper_map::const_iterator i = m->begin ();
    while (i != m->end ()) {
      i++;
      n++;
    }
    fi++;
  }
  vector <DOFMAPPER_BASE *> dofmappers_list(n); 
  sizet idofmappers = 0;
  fi = _dofmappers_map.begin ();
  while (fi != _dofmappers_map.end ()) {
    dofmapper_map *m = fi->second;
    dofmapper_map::const_iterator i = m->begin ();
    while (i != m->end ()) {
      DOFMAPPER_BASE *dofmapper = i->second;
      dofmappers_list[idofmappers] = dofmapper;
      idofmappers++; i++;
    }
    fi++;
  }
  CHECK(idofmappers == n, EXCEPTION_BAD_VALUE,;);
  sort (dofmappers_list.begin (), dofmappers_list.end (),
        EVS_OOOFS::DOFMAPPER_BASE_COMPARATOR());
  for (vector <DOFMAPPER_BASE *>::iterator i=dofmappers_list.begin();
       i != dofmappers_list.end(); i++) {
    DOFMAPPER_BASE *dofmapper = *i;
    const int neqns = dofmapper->neqns ();
    int *en = new int [neqns];
    bool *c = new bool [neqns];
    // first get the current equation constraints
    dofmapper->constrained (c);
    // now number the unconstrained equations: 1-based numbering!!!
    for (int j = 0; j < neqns; j++) {
      if (!c[j]) { en[j] = _total_of_equations+1; _total_of_equations++; }
      else       { en[j] = DOFMAPPER_BASE::invalid_eqnum; }
    }
    // and set the equation numbers in the dofmapper
    dofmapper->set_eqnums (en);
    delete [] en;
    delete [] c;
  }
  if (! _conn_graph_map.empty()) make_col_heights ();
}
#else
void EVS_OOOFS::naive_equation_numbering ()
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
      // now number the unconstrained equations: 1-based numbering!!!
      for (int j = 0; j < neqns; j++) {
        if (!c[j]) { en[j] = _total_of_equations+1; _total_of_equations++; }
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
  if (! _conn_graph_map.empty()) make_col_heights ();
}
#endif

int EVS_OOOFS::tot_of_equations ()
{
  if (_total_of_equations <= 0) {
    naive_equation_numbering ();
  }
  CHECK (_total_of_equations > 0, EXCEPTION_BAD_VALUE,;);
  return _total_of_equations;
}



bool EVS_OOOFS::solve (double shift)
{
  _eigvec.resize (_total_of_equations, _numberOfRequiredEigenValues); _eigvec.zero();
  _eigval.resize (_numberOfRequiredEigenValues); _eigval.zero();
  SMART_HANDLE<LOGGER_STREAM> lsb = _mgr->logger("solve", true);
  if (!_ssi.solve (*lsb, _K, _M,
                   &_eigval, &_eigvec, _rtolv, _max_iter,
		   _numberOfRequiredEigenValues)) {
    WARN("Not converged in SubspaceIteration");
  }
  *lsb << "done" << endl;
  return true;
}

void EVS_OOOFS::conn_graph_new_batch ()
{
  _making_conn_graph = true;
  _batch.clear();
}

void EVS_OOOFS::update_conn_graph ()
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

void EVS_OOOFS::add_to_batch (DOFMAPPER_BASE *dm)
{
  _batch.insert (dm);
}

void EVS_OOOFS::make_col_heights ()
{
  // Allocate and initialize the array
  _col_heights = new int [_total_of_equations];
  for (int k = 0; k < _total_of_equations; k++) _col_heights[k] = INT_MAX;

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
                _col_heights[eqnums[k]-1] = min(_col_heights[eqnums[k]-1], eqnums2[j]);
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

}


