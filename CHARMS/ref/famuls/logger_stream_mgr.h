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
#ifndef LOGGER_STREAM_MGR_H
#   define LOGGER_STREAM_MGR_H

#include<string>
#include "logger_stream.h"

class LOGGER_STREAM_MGR {

 public:

  LOGGER_STREAM_MGR () : _did_endl_last(true), _nblocks(0)
    {
      _to_stream = &std::cerr;
    }

  LOGGER_STREAM *logger (std::string where, bool print_exec_time);

  LOGGER_STREAM *logger (std::string where) {
    return this->logger (where, false);
  }

  void delete_notify (LOGGER_STREAM *logger) {
    _nblocks--;
  }
  
  void set_did_endl_last (bool did_endl_last) { _did_endl_last = did_endl_last; }
  bool did_endl_last () const { return _did_endl_last; }

 private: // private members ///////////////////////////
  
  bool          _did_endl_last;
  std::ostream *_to_stream;
  unsigned int  _nblocks;

};

#include "logger_stream_op.h"

#endif
