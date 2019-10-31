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
#include "bfun_set.h"
#include "field.h"
#include "algo_refine.h"
#include "ref_fen_src.h"
#include "ref_ctx.h"
#include "algo_err_rep.h"

ALGO_REFINE::ALGO_REFINE (string name, MGR *mgr,
                          GMESH *gmesh) : ALGO (name, mgr)
{
  string path = "algorithms/refine/" + this->name();
   
  _ref_fraction = 1.0 ;
  if (this->mgr()->db()->param_defined (path + "/ref_fraction")) {
    _ref_fraction = this->mgr()->db()->DB_GET_DOUBLE (path + "/ref_fraction");
  }
  
  _h_over_hbar_ref = 1;
  if (this->mgr()->db()->param_defined (path + "/h_over_hbar_ref")) {
    _h_over_hbar_ref = this->mgr()->db()->DB_GET_DOUBLE (path + "/h_over_hbar_ref");
  }
  //  CHECK(_h_over_hbar_ref >= 1, EXCEPTION_BAD_VALUE,;);
  
  _h_over_hbar_unref = 1; 
  if (this->mgr()->db()->param_defined (path + "/h_over_hbar_unref")) {
    _h_over_hbar_unref = this->mgr()->db()->DB_GET_DOUBLE (path + "/h_over_hbar_unref");
  }
  //  CHECK(_h_over_hbar_unref <= 1, EXCEPTION_BAD_VALUE,;); //RAW ashi
  
  _max_ref_level = 100000000; // unlimited
  if (this->mgr()->db()->param_defined (path + "/max_ref_level")) {
    _max_ref_level = this->mgr()->db()->DB_GET_INTEGER (path + "/max_ref_level");
  }

   _gmesh = gmesh;
  _true_hierarchical = false;
  if (this->mgr()->db()->param_defined (path + "/true_hierarchical")) {
    _true_hierarchical = this->mgr()->db()->DB_GET_BOOL (path + "/true_hierarchical");   
  } 
  _ref_ctx = new REF_CTX (_gmesh);
   
  _target_nbfuns = 0; // meaning *not specified*
  if (this->mgr()->db()->param_defined (path + "/target_nbfuns")) {
    _target_nbfuns = this->mgr()->db()->DB_GET_INTEGER (path + "/target_nbfuns");   
  } 
}

ALGO_REFINE::~ALGO_REFINE()
{
  if (_ref_ctx) delete _ref_ctx;
  
}


bool ALGO_REFINE::check_rule(int rule,BFUN *bfun)
{
  switch (rule)
  {
    case 1:
       return 1;
    case 2:
//rule for refinement
      return (rule_2(bfun));
    case 3:
// rule for unrefinement
      return (rule_3(bfun));
    default:
      return (0);
  }
  return (1);
}


bool ALGO_REFINE::rule_2 (BFUN *bfun)
{
  // a function may be refined only when all its parents are refined
  set<FEN*> parent_set;
  bfun->unrefinement_set(&parent_set);
  for (set<FEN*>::iterator i = parent_set.begin(); i != parent_set.end(); i++) {
    BFUN *bfunx = bfun->bfun_set()->fen_to_bfun_any(*i);
    if (!bfunx->is_refined()) return false;
  }
  return true;
}

bool ALGO_REFINE::rule_3 (BFUN *bfun)
{
// A function on level j may be unrefined  only if 
            // a. it was previously refined  
            // b. all its children on level j+1 are not refined
  if(!bfun->is_refined()) return false;
  set <FEN*> child_set;
  bfun->refinement_set(&child_set);  
  for (set <FEN*>::iterator i= child_set.begin(); i != child_set.end(); i++) {
    BFUN *bfun1 = bfun->bfun_set()->fen_to_bfun_any(*i);
    if (bfun1 != 0)
      if (bfun1->is_refined()) return false;
  }   
  return true;
}






