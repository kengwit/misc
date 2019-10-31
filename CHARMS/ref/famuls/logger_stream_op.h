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

#ifndef LOGGER_STREAM_OP_H
# define LOGGER_STREAM_OP_H 1

template <class T>
inline
LOGGER_STREAM &
LOGGER_STREAM::operator<< (const T& t)
{
  if (_logger_stream_mgr->did_endl_last()) {
    do_header ();
    _logger_stream_mgr->set_did_endl_last (false);
  }

  *_to_stream << t;

  return *this;
}

#endif
