#include "watchpoint_ecell.h"
#include "timeP.h"
#include <fstream>
#include <ctime>
#include "ecell.h"


bool WATCH_POINT_ECELL::set_ecell (ECELL *ecell) {
  _ecell = ecell;
  if (ecell) { //unmap_watchpoint sets ecell to 0 
    GCELL *gcell = ecell->gcell()->root();
    POINT param_loc = _input_param_loc;
    while (gcell  != ecell->gcell() ) {    
      GCELL *child_gcell;
      if (!gcell->map_to_child (param_loc,&child_gcell,&_param_loc)) return false;
      gcell = child_gcell;
      param_loc = _param_loc;
      CHECK(gcell, EXCEPTION_NULL_PTR,;);
    }
    _param_loc = param_loc;
  }
  return true;
} 

ECELL *WATCH_POINT_ECELL::ecell() {return _ecell;}


