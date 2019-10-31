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
#include "gcell.h"
#include "field.h"

OFACTORY<GCELL, GCELL::make_gcell_func, vector <FEN *> > *GCELL::_gcell_factory = 0;

GCELL::GCELL ()
{
  _gcell_group = 0;
  _id          = ID_NONE;
}

bool
GCELL::register_make_func (make_gcell_func make, string gcell_type, string implementation)
{
  if (_gcell_factory == 0) {
    _gcell_factory = new OFACTORY<GCELL, make_gcell_func, vector <FEN *> >;
  }
  return (GCELL::_gcell_factory)->register_make_func (make, gcell_type, implementation);
}


GCELL *
GCELL::make_gcell (string type, string implementation, vector <FEN *> fens)
{
  CHECK (_gcell_factory, EXCEPTION_NULL_PTR,;);
  return (GCELL::_gcell_factory)->make (fens, type, implementation);
}


GCELL *
GCELL::map_to_topmost_gcell_in_field (FIELD_BASE *field,
                                      POINT &param_loc, POINT *mapped_to_param_loc)
{
  POINT curr_param_loc = param_loc;
  GCELL *child;
  GCELL *tested = this;        // start with self
  // Try to get to the top of the pyramid
  POINT child_param_loc;
  while (tested->map_to_child (curr_param_loc, &child, &child_param_loc)) {
    tested = child;
    curr_param_loc = child_param_loc;
  }
  // Now we should be at the top, descend.  Once we find
  // a cell used by active basis functions in the field, we use it.
  while (! field->active_over (tested)) {
    GCELL *parent;
    POINT parent_param_loc;
    if (!tested->map_to_parent (curr_param_loc, &parent, &parent_param_loc)) {
      return 0; // can't go deeper, and still haven't found a usable cell
    }
    curr_param_loc = parent_param_loc;
    tested = parent;
  }
  // We've got it.
  *mapped_to_param_loc = curr_param_loc;
  return tested;
}

bool
GCELL::is_topmost_gcell_in_field (FIELD_BASE *field)
{
  // Is the field active over self?  If not, there is no point in testing further.
  if (!field->active_over (this)) return false;
  // Does this cell have any children?  If not, there is no point in further tests.
  if (this->nchildren () == 0) return true;
  
  // Now test if there is a descendant of this gcell which is active in field.
  // If so, this is not the topmost cell.
  stack <GCELL *> s; 
  s.push (this);
  while (! s.empty ()) {
    GCELL *tested = s.top(); s.pop();
    for (sizet j = 0; j < tested->nchildren (); j++) {
      GCELL *child = tested->child(j);
      if (field->active_over (child)) return false; // descendant active: done
      s.push (child);
    }
  }
  return true; // no descendant active: done
}

void
GCELL::debug_display ()
{
  cerr << "GCELL";
}

sizet
GCELL::level ()
{
  sizet l = 0;
  GCELL *gcell = this;
  while (gcell->parent()) {
    l++;
    gcell = gcell->parent();
  }
  return l;
} 

double
GCELL::char_dim ()
{
  const sizet nfens = conn()->nfens();
  double maxdist2 = 0;
  for (sizet l=0; l<nfens; l++) {
    POINT a = conn()->fen(l)->ref_loc();
    for (sizet k=0; k<nfens; k++) {
      POINT b = conn()->fen(k)->ref_loc();
      double dist2 = ((a[0]-b[0])*(a[0]-b[0]) +(a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]));
      if (dist2 > maxdist2) maxdist2 = dist2;
    }
  }
  return sqrt(maxdist2); 
}
