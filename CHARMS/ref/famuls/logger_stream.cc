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
#include "logger_stream.h"
#include "logger_stream_mgr.h"

LOGGER_STREAM::~LOGGER_STREAM() {
  _logger_stream_mgr->delete_notify (this);
}

LOGGER_STREAM *
LOGGER_STREAM::clone (std::string where) {
  return _logger_stream_mgr->logger(where, _print_exec_time);
}

void
endl (LOGGER_STREAM &s)
{
  if (s._print_exec_time) {
    double ut = s._timer.lap();
    *s._to_stream << " [" << ut << " sec]";
  }
  *s._to_stream << std::endl;
  s._logger_stream_mgr->set_did_endl_last (true);
}

void
endl_notime (LOGGER_STREAM &s)
{
  *s._to_stream << std::endl;
  s._logger_stream_mgr->set_did_endl_last (true);
}

