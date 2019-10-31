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
#ifndef LOGGER_STREAM_H
#   define LOGGER_STREAM_H

#include <string>
#include <iostream>
#include "timer.h"

class LOGGER_STREAM_MGR;

/**
   Log information during program execution.
 */
class LOGGER_STREAM {
  
 public: // object methods ///////////////////////////////////
  
  /**
   */
  LOGGER_STREAM (class LOGGER_STREAM_MGR *logger_stream_mgr,
                 std::ostream *to_stream, unsigned int indent, bool print_exec_time) :
    _logger_stream_mgr(logger_stream_mgr),
    _to_stream(to_stream),
    _indent(indent),
    _print_exec_time(print_exec_time)
    {
    }

  LOGGER_STREAM *clone (std::string where);
  
  void set_indent (unsigned int indent) { _indent = indent; }

  ~LOGGER_STREAM ();

  template <typename T>
    LOGGER_STREAM &operator<< (const T &t);

  LOGGER_STREAM &operator<< (void (func)(LOGGER_STREAM &s)) {
    func (*this);
    return *this;
  }
  
  friend void endl (LOGGER_STREAM &s);
  friend void endl_notime (LOGGER_STREAM &s);
  
 private: // object data ////////////////////////////////////////
  
  class LOGGER_STREAM_MGR   *_logger_stream_mgr;
  std::ostream              *_to_stream;
  unsigned int               _indent;
  bool                       _print_exec_time;
  TIMER                      _timer;
  
 private: // object methods /////////////////////////////////////
  
  void do_header () {  
    *_to_stream << "#";
    for (unsigned int j = 0; j < _indent; j++) {
      *_to_stream << "   ";
    }
  }
  
};

extern void
endl (LOGGER_STREAM &s);

extern void
endl_notime (LOGGER_STREAM &s);


#endif

