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
#ifndef ELIXIR_INTERFACE_H
# define ELIXIR_INTERFACE_H

#include "famuls.h"
#include "famexception.h"

#if defined (GIFACE) && GIFACE
extern "C" {
#include "Esimple.h"
}
#endif

class ELIXIR_INTERFACE {

 public:
  
  static const sizet MAX_COLORS = 13;
  
  static const char *color_names_table[MAX_COLORS]; 

  public:

#if defined (GIFACE) && GIFACE
  static void build_interface (void *client, void (*callback)(void *client, Widget parent_pane));
#endif

 private:

  ELIXIR_INTERFACE () {}

  static bool initialized;

 public:

#if defined (GIFACE) && GIFACE
  static void *client;
  static void (*callback)(void *client, Widget parent_pane);
#endif

};

#endif

