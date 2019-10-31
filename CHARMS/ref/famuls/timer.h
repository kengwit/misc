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
#ifndef TIMER_H
# define TIMER_H

#include "timeP.h"

class TIMER {
  
 public:
  
  TIMER () {
    _timer = ckit_timer_new ();
  }
  void start () {
    ckit_timer_restart (_timer);
  }    
  double lap () {
    double ut = ckit_timer_lap (_timer);
    this->start ();
    return ut;
  }    
  ~TIMER () {
    ckit_timer_free (_timer);
  }
  
 private:
  
  ckit_timer_t *_timer;

};

#endif
