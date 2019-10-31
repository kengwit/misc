#include "field.h"

template <int NUM_COMPONENTS>
size_t
NZ::nzs (BFUN_DOFPARAM_PAIR_ID dpid)
{
  set <BFUN_DOFPARAM_PAIR_ID> ninteractionset;
  ninteractionset.clear(); ninteractionset.insert (dpid); // self
  FIELD_PAIR<NUM_COMPONENTS> *fp = _f->field_pair (dpid);
  BFUN *bfun = fp->bfun();
  size_t ngcells = bfun->ngcells();
  for (size_t j = 0; j < ngcells; j++) {
    GCELL *gcell = bfun->gcell(j);
    set <BFUN_DOFPARAM_PAIR_ID> actpairs = _f->pairs_active_over_gcell (gcell);
  }
  return _nzs;
}
