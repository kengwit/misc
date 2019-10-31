#ifndef WATCH_POINT_ECELL_H
# define WATCH_POINT_ECELL_H

#include "gcell.h"
#include "watchpoint.h"

class WATCH_POINT_ECELL : public WATCH_POINT {
  
 public: 
  
  WATCH_POINT_ECELL (list<string> var, string prefix, sizet gcell_id, POINT input_param_loc)
    : WATCH_POINT(var, prefix) {
    _fen_id          = -1;
    _gcell_id        = gcell_id;
    _input_param_loc = input_param_loc;
    _param_loc       = input_param_loc;
    _ecell           = 0;
  }

  ~WATCH_POINT_ECELL () {
  }
  
  POINT  param_loc() {return _param_loc;}
  sizet gcell_id() {return _gcell_id;} 
  void  map_to_ecell (class ECELL *ecell, POINT param_loc);
  bool set_ecell (class ECELL *ecell) ; 
  class ECELL *ecell() ;
  
 private:
  
  sizet          _gcell_id;
  long int       _fen_id;
  POINT          _input_param_loc;
  class ECELL   *_ecell;
  POINT          _param_loc;
  
};

#endif
