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
#include "bfun_set.h"
#include "gcell.h"
#include "gmesh.h"
#include "gsubmesh.h"
#include "field.h"
#include "evalpt.h"
#include <map>
#include "proto_field_transfer.h"



void
BFUN_SET::build_from_gcells (list <GCELL *> igcells,  BFUN_SET *source_bfun_set)
{
  vector <BFUN*> inactive_bfun; inactive_bfun.clear();
  typedef   map <FEN *, vector<GCELL*> > fen_to_gcells_map_t;
  // build up a map from nodes to gcells sharing those nodes
  fen_to_gcells_map_t fen_to_gcells_map;
  for (list <GCELL *>::const_iterator gli = igcells.begin (); gli != igcells.end (); gli++) {
    GCELL *gcell = *gli;
    const CONN_BASE *conn = gcell->conn ();
    sizet nfens = conn->nfens ();
    for (unsigned int j = 0; j < nfens; j++) {
      FEN *fen = conn->fen(j);
      fen_to_gcells_map_t::iterator i = fen_to_gcells_map.find (fen);
      if (i == fen_to_gcells_map.end ()) {
        vector<GCELL*> gv;
        fen_to_gcells_map.insert (fen_to_gcells_map_t::value_type (fen, gv));
        i = fen_to_gcells_map.find (fen);
      }
      (i->second).push_back (gcell);
    }
  }
  // now iterate the map and set up the basis functions
  
  for (fen_to_gcells_map_t::const_iterator i = fen_to_gcells_map.begin ();
       i != fen_to_gcells_map.end (); i++) {
    FEN *fen = i->first;
    // Decide whether that node is to be included:
    vector <GCELL *> gcells = i->second;
    BFUN_FE *bfun = new BFUN_FE (fen, gcells,this); //RAW : should be replaced with a make function.
    bfun->set_is_active(true);
    this->add (bfun);
  }     
}


void BFUN_SET::debug_display(char filename[])
{
}

void BFUN_SET::debug_display()
{
  cerr << "BFUN_SET" << endl;
  cerr <<"_nbfuns_total = "<< _nbfuns_total << "\n";
  cerr <<"_nbfuns_active = "<< _nbfuns_active << "\n";
  for (sizet i = 0; i < _nbfuns_total; i++) {
    CHECK(_bfuns[i]->bfun_set() == this, EXCEPTION_BAD_VALUE,;);
    cerr << "fen " << _bfuns[i]->fen()->id() << (_bfuns[i]->is_active() ? " (active)" : "") << endl;   
  }   
}


list <GCELL *>
BFUN_SET::build_gcell_list ()
{
  list <GCELL *> igcells;
  list <GSUBMESH *>::const_iterator si = _gsubmeshes.begin ();
  while (si != _gsubmeshes.end ()) {
    GSUBMESH *gsubmesh = *si;
    GSUBMESH::gcell_group_enumerator_t e = gsubmesh->gcell_group_enumerator ();
    GCELL_GROUP *gcell_group;
    e.reset ();
    while ((gcell_group = e.next ())) {
      GCELL_GROUP::gcell_enumerator_t e1 = gcell_group->gcell_enumerator ();
      GCELL *gcell;
      e1.reset ();
      while ((gcell = e1.next ())) {
        igcells.push_back (gcell);
      }
    }
    si++;
  }
  return igcells;
}

BFUN_SET::BFUN_SET (REFINEMENT_STRATEGY refinement_strategy )
{
  
  _in_constructor = true;
  _refinement_strategy = refinement_strategy;
  _gmesh = 0;
  _gsubmeshes.clear ();
  _bfuns.clear();
  _nbfuns_active = _nbfuns_total = 0;
  _bfun_dofparam_ids.clear ();
  _in_constructor = false;
  _mem_raise_size = 10;
#if defined(USE_PU) && USE_PU
  _is_pu = false;
#endif
}

/* set of basis functions on the domain defined by a set of gsubmeshes */
BFUN_SET::BFUN_SET (REFINEMENT_STRATEGY refinement_strategy,list <GSUBMESH *> gsubmeshes)
{
  _in_constructor = true;  
  _refinement_strategy = refinement_strategy;
  _gsubmeshes.clear ();
  _gsubmeshes.insert (_gsubmeshes.begin (), gsubmeshes.begin (), gsubmeshes.end ());
  GSUBMESH *first_submesh = *_gsubmeshes.begin();
  CHECK (first_submesh, EXCEPTION_NULL_PTR,;);
  _gmesh = first_submesh->gmesh(); // `local' field: constructed from submeshes
  _bfuns.clear();
  _nbfuns_active = _nbfuns_total = 0; 
  _bfun_dofparam_ids.clear ();
  
  list <GCELL *> igcells = build_gcell_list ();
  CHECK (!igcells.empty (), EXCEPTION_BAD_VALUE,;);
  
  build_from_gcells (igcells, 0);
  _mem_raise_size = 10;
#if defined(USE_PU) && USE_PU
  _is_pu = false;
#endif
  _in_constructor = false;
}

/* set of basis functions on the domain defined by a set of gsubmeshes */
BFUN_SET::BFUN_SET (REFINEMENT_STRATEGY refinement_strategy, GMESH *gmesh)
{
  _in_constructor = true;
  _refinement_strategy = refinement_strategy;  
  _gmesh = gmesh; // `global' field: constructed from the mesh
  _gsubmeshes.clear ();
  GMESH::gsubmesh_enumerator_t e = gmesh->gsubmesh_enumerator ();
  GSUBMESH *gsubmesh;
  e.reset ();
  while ((gsubmesh = e.next ())) {
    _gsubmeshes.push_back (gsubmesh);
  }
  _bfuns.clear();
  _nbfuns_active = _nbfuns_total = 0;
  _bfun_dofparam_ids.clear ();
  
  list <GCELL *> igcells = build_gcell_list ();
  CHECK (!igcells.empty (), EXCEPTION_BAD_VALUE,;);

  build_from_gcells (igcells, 0);
  _mem_raise_size = 10;
#if defined(USE_PU) && USE_PU
  _is_pu = false;
#endif
  _in_constructor = false;
}


BFUN *
BFUN_SET::fen_to_bfun_any (FEN *fen)
{
  BFUN_DOFPARAM_PAIR_ID dpid = this->dofparam_id (fen);
  if (dpid != INVALID_BFUN_DOFPARAM_PAIR_ID){
     BFUN *bfun = _bfuns[dpid];
     CHECK(bfun->fen()->id()==fen->id(),EXCEPTION_BAD_VALUE,;);//RAW
     if (bfun->bfun_set() == this) return bfun;
     else                          return 0;
  } 
  else  return 0;
}


BFUN *
BFUN_SET::fen_to_bfun (FEN *fen)
{
  BFUN* bfun = fen_to_bfun_any (fen);
  if (bfun != 0) {   
    if (bfun->is_active()) return bfun; else return 0;
  } 
  else return 0;
}



void
BFUN_SET::display ()
{
  {
    cerr << "Gsubmeshes: " << endl;
    GSUBMESH *gsubmesh;
    BFUN_SET::gsubmesh_enumerator_t e = gsubmesh_enumerator ();
    e.reset ();
    while ((gsubmesh = e.next ())) {
      cerr << gsubmesh->name () << endl;
    }
  }
  {
    cerr << "Bfuns: " << endl;
    BFUN *bfun;
    BFUN_SET::bfun_enumerator_t e = bfun_enumerator ();
    e.reset ();
    while ((bfun = e.next ())) {
      cerr << "Node " << bfun->fen()->id() << endl;
    }
  }
}

SMART_HANDLE <BFUN_SET> BFUN_SET::clone ()
{
  SMART_HANDLE <BFUN_SET> new_bfun_set = new BFUN_SET (this->_refinement_strategy);
  new_bfun_set->_in_constructor = true; // paired with commit clone.
  new_bfun_set->_gmesh = this->_gmesh;
  new_bfun_set->_gsubmeshes.resize(_gsubmeshes.size()); //RAW
  copy (_gsubmeshes.begin (), _gsubmeshes.end (), new_bfun_set->_gsubmeshes.begin ());
 // cout <<"list of added bfuns\n"; 
  for (sizet bfun_count = 0; bfun_count<_nbfuns_total; bfun_count++) {
    BFUN *new_bfun = _bfuns[bfun_count]->clone(new_bfun_set.get_ptr());
    new_bfun->set_is_active(_bfuns[bfun_count]->is_active());
    new_bfun_set->add (new_bfun);
   // cout<<"added bfun"<<new_bfun->fen()->id()<<"\n"; /
  }
  new_bfun_set->_bfun_dofparam_ids.resize (this->_bfun_dofparam_ids.size ());
  copy (this->_bfun_dofparam_ids.begin (), this->_bfun_dofparam_ids.end (),
        new_bfun_set->_bfun_dofparam_ids.begin ());
#if defined(USE_PU) && USE_PU
  new_bfun_set->_is_pu = _is_pu;
#endif
  return new_bfun_set;
}

bool BFUN_SET::is_present(BFUN*bfun) {
  return (this->dofparam_id(bfun->fen()) != INVALID_BFUN_DOFPARAM_PAIR_ID) ; 
} 

void BFUN_SET::swap(sizet temp_loc1,sizet temp_loc2)
{
  BFUN *temp_bfun   = _bfuns[temp_loc1];
  _bfuns[temp_loc1] = _bfuns[temp_loc2];
  _bfuns[temp_loc2] = temp_bfun;
  
  BFUN_DOFPARAM_PAIR_ID temp_dpid = _bfun_dofparam_ids[_bfuns[temp_loc1]->fen()->uniqobjid()];
  _bfun_dofparam_ids[_bfuns[temp_loc1 ]->fen()->uniqobjid()]
    = _bfun_dofparam_ids[_bfuns[temp_loc2]->fen()->uniqobjid()]; 
  _bfun_dofparam_ids[_bfuns[temp_loc2]->fen()->uniqobjid()] = temp_dpid;
}

void BFUN_SET::insert_bfun (BFUN *bfun)
{
  _mem_raise_size = 10;
  if (_bfuns.size()<=_nbfuns_total) {
    _bfuns.resize(_bfuns.size()+_mem_raise_size);   
  }  
  _bfuns[_nbfuns_total] = bfun;
  _nbfuns_total++;

  if(_bfun_dofparam_ids.size()<= bfun->fen()->uniqobjid()) {
     sizet p_size = _bfun_dofparam_ids.size();
    _bfun_dofparam_ids.resize(bfun->fen()->uniqobjid()+_mem_raise_size);
    for (sizet k = p_size; k< _bfun_dofparam_ids.size(); k++) {
      _bfun_dofparam_ids[k] = INVALID_BFUN_DOFPARAM_PAIR_ID;   
    }
  }  
  _bfun_dofparam_ids[bfun->fen()->uniqobjid()] = _nbfuns_total-1;   
}

void
BFUN_SET::add (BFUN *bfun)
{   
  CHECK (_in_constructor,EXCEPTION_ILLEGAL_USE,;); //RAW : now we are adding bfun in clone 
  CHECK(bfun->bfun_set()== this, EXCEPTION_ILLEGAL_USE,;);
  if (!is_present(bfun)) insert_bfun(bfun);
  if (bfun->is_active()) { //if basis function has to be  active. 
    if(!in_active_set(bfun)){ // but it was grouped among inactive
      swap(dofparam_id(bfun->fen()),_nbfuns_active);
      _nbfuns_active++;
    } 
  } else { // if a basis function has to be inactive 
    if (in_active_set(bfun)) { //but it was grouped among actives 
      swap(dofparam_id(bfun->fen()), _nbfuns_active-1); 
      _nbfuns_active--;
    }
  } 
}

 
BFUN_DOFPARAM_PAIR_ID
BFUN_SET::dofparam_id (FEN *fen)
{
  if (fen->uniqobjid() >= _bfun_dofparam_ids.size ())
    return INVALID_BFUN_DOFPARAM_PAIR_ID;
  else
    return _bfun_dofparam_ids[fen->uniqobjid()];
}


BFUN_SET::~BFUN_SET () {
  // RAW delete the basis functions
  for (vector <BFUN *>::iterator i = _bfuns.begin();
       i != _bfuns.end(); i++) {
    BFUN *bfun = *i;
    delete bfun;
  }
}

BFUN_DOFPARAM_PAIR_ID
BFUN_SET::dofparam_id (BFUN *bfun)
{
  CHECK (bfun->bfun_set() == this, EXCEPTION_BAD_ACCESS,;);
  if (!bfun->is_active()) return  INVALID_BFUN_DOFPARAM_PAIR_ID;//RAW
  else                    return  dofparam_id (bfun->fen());
}

void BFUN_SET::commit_clone()
{
#if defined(USE_PU) && USE_PU
  make_pu();
#endif
  _in_constructor = false;
}


#if defined(USE_PU) && USE_PU
void BFUN_SET::make_pu() 
{
  for (vector <BFUN *>::iterator i = _bfuns.begin(); i != _bfuns.end(); i++) {
    BFUN *bfun = *i;
    if (bfun) bfun->set_beta(1);
  } 
  _is_pu = false;
  if (_refinement_strategy == QUASI_HIERARCHICAL) {
    SMART_HANDLE<BFUN_SET> bfun_set(this);
    FIELD_SCALAR *betas = new FIELD_SCALAR("betas", bfun_set, 0);
    PROTO_FIELD_TRANSFER <1> p;
 //  cerr << "before  transfer @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl; //RAW ashi
 //   betas->debug_display(); //RAW ashi

    p.transfer(1,betas);
  //  p.check_pu (betas, 47, 1000 ); 
 //   cerr << "after transfer @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl; //RAW:ashi 
//    betas->debug_display(); //RAW :ashi
    for (sizet  i = 0; i <  betas->npairs(); i++) {
      FIELD_PAIR <1> *fp  = betas->field_pair(i);
/*
  In case of quadratic elements some of the basis functions gets omitted
despite following quadratic consistency by the nature of mature refinement algorithm. Following block deactivates those unnessasry basis functions.
*/
// cout << "following basis functions are forced to be inactive \n" ; //RAW
      if (fp->dofparam()(0) != 0) {
        fp->bfun()->set_beta(fp->dofparam()(0));
      } else {
        fp->bfun()->set_is_active(false);
        this->add(fp->bfun());
        cout << fp->bfun()->fen()->id() << "\t is forced to be inactive \n" ; //RAW  
      }    
    } 
    delete betas; 
    bfun_set.forget_me(this);// prevent the handle from deleting `this'
    _is_pu = true;
  }
  //this->debug_display("in make_pu()");
}
/*
void  BFUN_SET::check_consist  (FIELD<NUM_COMPONENTS> *geometry)
{
      PROTO_FIELD_TRANSFER <1> p;
      sizet maxlevel = p.max_level (geometry);
      for (sizet level = 0; level < maxlevel+1; level++) {
        const sizet npairs = geometry->npairs();
        for (sizet j = 0; j < npairs; j++) {
          FIELD_PAIR<NUM_COMPONENTS> *fp = geometry->jth_field_pair (j);
          BFUN *bfun = fp->bfun();
          if (bfun->level() == level) {
            FEN *fen = bfun->fen();
            GCELL *gcell = bfun->gcell(0); // any gcell would do
            EVALPT geomevalpt (geometry, gcell, fen);
            geometry.eval ();
            FIXED_VECTOR<NUM_COMPONENTS> srcv(sval);
            sizet nbfuns = geomevalpt.nbfuns();
            for (sizet J = 0; J < nbfuns; J++) {
                double N = geomevalpt.N(J);
                BFUN_DOFPARAM_PAIR_ID geomdpid = geometry->dofparam_id (bfun->fen());
                CHECK (destdpid != INVALID_BFUN_DOFPARAM_PAIR_ID, EXCEPTION_NULL_PTR,;);
                CHECK (geomdpid != INVALID_BFUN_DOFPARAM_PAIR_ID, EXCEPTION_NULL_PTR,;);
                //calculate summation (n*beta*x), it should be equal to x for linear elements and x^2 for quadratic elements 
                FIXED_VECTOR<NUM_COMPONENTS> dp = bfun->beta() * geometry->field_pair(geomdpid)->dofparam();
                srcv.add (N, dp);
            }
            CHECK(N_J != 0, EXCEPTION_BAD_VALUE,;); // RAW: expensive, use only for debug
            srcv.scale(1/N_J); // Adjust for PU
            fp->set_dofparam_all (srcv);
          } // if
        } // for
      } // for
}
*/

#endif


