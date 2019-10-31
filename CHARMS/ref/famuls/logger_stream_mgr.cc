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

#include "logger_stream_mgr.h"

LOGGER_STREAM *
LOGGER_STREAM_MGR::logger (std::string where, bool print_exec_time) {
  LOGGER_STREAM *ls = new LOGGER_STREAM (this, _to_stream, _nblocks++, print_exec_time);
  if (!_did_endl_last) *ls << endl;
  *ls << where << endl_notime;
  ls->set_indent (_nblocks);
  return ls;
}
