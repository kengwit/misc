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
#ifndef ALGO_INTERPOLATE_H
# define ALGO_INTERPOLATE_H

#include "algo.h"
#include "field.h"
#include "les.h"
#include "gcell.h"
#include "evalpt.h"

template <int NUM_COMPONENTS>
class ALGO_INTERPOLATE: public ALGO {

 public:

  ALGO_INTERPOLATE ();
  ~ALGO_INTERPOLATE () {}

  void interpolate (FIELD<NUM_COMPONENTS> *source,
                    FIELD_VECTOR *source_geometry,
                    FIELD<NUM_COMPONENTS> *destination,
                    FIELD_VECTOR *destination_geometry);

};

template <int NUM_COMPONENTS>
ALGO_INTERPOLATE<NUM_COMPONENTS>::ALGO_INTERPOLATE() : ALGO("") {
}

template <int NUM_COMPONENTS>
void ALGO_INTERPOLATE<NUM_COMPONENTS>::interpolate (FIELD<NUM_COMPONENTS> *source,
                                                    FIELD_VECTOR *source_geometry,
                                                    FIELD<NUM_COMPONENTS> *destination,
                                                    FIELD_VECTOR *destination_geometry
                                                    )
{
  char *argva[] = {
    "?", "-ksp_type", "cg", "-pc_type", "jacobi", 0
      };
  char **argv = argva;
  int argc = sizeof (argva) / sizeof (argva[0]) - 1;
  const string field_name = "interpolate_aux_field";
  FIELD_SCALAR aux (field_name, destination->bfun_set (), 0);
  LES les (string("les"), this->mgr(), string("fles_solver dmeschach"));

  // Set up the dofmappers
  typedef map <BFUN *, DOFMAPPER<1> *> dmmap_t;
  dmmap_t dmmap;
  {
    BFUN_SET *bfun_set = aux.bfun_set ();
    BFUN_SET::bfun_enumerator_t e = bfun_set->bfun_enumerator ();
    BFUN *bfun;
    e.reset ();
    while (bfun = e.next ()) {
      BFUN_DOFPARAM_PAIR_ID dpid = aux.dofparam_id (bfun->fen());
      DOFMAPPER<1> *dm = les.dofmapper (field_name, &aux, dpid);
      dmmap.insert (dmmap_t::value_type (bfun, dm));
    }
  }
  
  // Assemble the lhs
  les.zero_lhs ();
  les.dump_lhs ("lhs");
  {
    BFUN_SET *bfun_set = aux.bfun_set ();
    BFUN_SET::bfun_enumerator_t e = bfun_set->bfun_enumerator ();
    BFUN *bfun;
    e.reset ();
    while (bfun = e.next ()) {  // For each basis function
      // cerr << "Evaluating at " << bfun->fen()->id() << endl;
      // Evaluate the bfun set of the destination field (it shares the bfun_set with aux) 
      // at this fen: bfun->fen()
      CHECK (bfun->ngcells() >= 1, EXCEPTION_BAD_ACCESS,;);
      GCELL *gcell = bfun->gcell(0);
      POINT param_loc;
      gcell->map_fen (bfun->fen(), &param_loc); // map the fen ...
      double mapped_to_param_loc; // ... and find the topmost gcell in field containing this point
      GCELL *evgcell = gcell->find_topmost_gcell_in_field (&aux, param_loc, &mapped_to_param_loc);
      CHECK (evgcell, EXCEPTION_BAD_VALUE,;);
      EVALPT evalpt (&aux, destination_geometry, gcell, mapped_to_param_loc);
      evalpt.eval ();
      // 
      DOFMAPPER<1> *dofmap;
      {dmmap_t::const_iterator i = dmmap.find (bfun); CHECK (i != dmmap.end (), EXCEPTION_BAD_VALUE,;); dofmap = i->second;}
      // Assemble lhs
      double values[1][1];
      for (int k = 0; k < evalpt.nbfuns (); k++) {
        BFUN *obfun = evalpt.bfun (k);
        DOFMAPPER<1> *odofmap;
        {dmmap_t::const_iterator i = dmmap.find (obfun); CHECK (i != dmmap.end (), EXCEPTION_BAD_VALUE,;); odofmap = i->second;}
        values[0][0] = evalpt.N(k);
        // cerr << "   Interaction with " << obfun->fen()->id() << " (" << values[0][0] << ") " << endl;
        {
          int eqnums[1];
          dofmap->eqnums (eqnums);
          int oeqnums[1];
          odofmap->eqnums (oeqnums);
          // cerr << "     Adding " << values[0][0] << " in row " << eqnums[0] << " col " << oeqnums[0] << endl;
        }
        les.lhs_add (dofmap, odofmap, values);
        les.dump_lhs ("lhs");
      }
    }
  }
  les.finish_lhs ();
  les.dump_lhs ("lhs");

  // Then for each rhs
  for (sizet k = 0; k < NUM_COMPONENTS; k++) {
    les.zero_rhs ();
    BFUN_SET *bfun_set = aux.bfun_set ();
    BFUN_SET::bfun_enumerator_t e = bfun_set->bfun_enumerator ();
    BFUN *bfun;
    e.reset ();
    while (bfun = e.next ()) {
      // Evaluate the source field
      CHECK (bfun->ngcells() >= 1, EXCEPTION_BAD_ACCESS,;);
      GCELL *gcell = bfun->gcell(0);
      FIXED_VECTOR<NUM_COMPONENTS> v = source->evaluate (source_geometry, bfun->fen(), gcell);
      // Evaluate the bfun set of the destination field
      double xi, eta, theta;
      gcell->map_fen (bfun->fen(), &xi, &eta, &theta); // map the fen ...
      double xi_mapped, eta_mapped, theta_mapped; // ... and find the topmost gcell in field containing this point
      GCELL *evgcell = gcell->find_topmost_gcell_in_field (&aux, xi, eta, theta,
                                                           &xi_mapped, &eta_mapped, &theta_mapped);
      CHECK (evgcell, EXCEPTION_BAD_VALUE,;);
      EVALPT evalpt (&aux, destination_geometry, gcell, xi_mapped, eta_mapped, theta_mapped);
      evalpt.eval ();
      // 
      DOFMAPPER<1> *dofmap;
      {dmmap_t::const_iterator i = dmmap.find (bfun); CHECK (i != dmmap.end (), EXCEPTION_BAD_VALUE,;); dofmap = i->second;}
      double values[1];
      values[0] = v(k);
      les.rhs_add (dofmap, values);
    }
    les.finish_rhs ();
    les.dump_rhs ("rhs");
    les.solve ();
    les.solution (field_name, &aux);
    // Transfer aux field to the destination field
    FIELD_PAIR<NUM_COMPONENTS> *dp;
    FIELD_PAIR<1> *ap;
    for (sizet j = 0; j < destination->npairs (); j++) {
      ap = aux.field_pair (j);
      dp = destination->field_pair (j);
      CHECK (dp->bfun()->fen() == ap->bfun()->fen(), EXCEPTION_BAD_VALUE,;); // RAW not really necessary: they both use the same bfun_set
      // cerr << "Setting " << dp->bfun()->fen()->id() << " " << ap->dofparam_comp (0) << endl;
      dp->set_dofparam_comp (k, ap->dofparam_comp (0));
    }
  }
}


#endif
