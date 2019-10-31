/*
                           OOOFS
                 Copyright (C) 2003, Petr Krysl
                 
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/
#ifndef ooofs_error_h
# define ooofs_error_h

#include <iostream>

namespace OOOFS {
  typedef int error_notify_func (int error_num,
                                 int line,
                                 char *file,
                                 char *detail,
                                 void *client_data);
  extern error_notify_func *error_notify;

#undef _error
#define _error(s)                                       \
{                                                       \
  if (OOOFS::error_notify)                              \
    OOOFS::error_notify(0, __LINE__, __FILE__, s, 0);   \
  else {                                                \
    std::cerr << "OOOFS error: " << s << std::endl;               \
  }                                                     \
}
  extern error_notify_func *error_notify;
}
#endif
