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
#include "elixir_interface.h"
#include "args.h"
#if defined (GIFACE) && GIFACE
bool ELIXIR_INTERFACE::initialized = false;
void *ELIXIR_INTERFACE::client = 0;
void (*ELIXIR_INTERFACE::callback)(void *client, Widget parent_pane) = 0;

void
ELIXIR_INTERFACE::build_interface (void *client, void (*callback)(void *client, Widget parent_pane)) {
  if (!ELIXIR_INTERFACE::initialized) {
    ELIXIR_INTERFACE::client = client;
    ELIXIR_INTERFACE::callback =  callback;
    unsigned long m = ESI_GRAPHIC_EDITOR_MASK;
    ESIBuildInterface (m, ARGS::argc, ARGS::argv);
    EMSetClient (ESIModel (), client);
    ESIHandleCmd ("color red");
    ESIHandleCmd ("shrink 1.0");
    ESIHandleCmd ("fill solid");
    ELIXIR_INTERFACE::initialized = true;
  }
}

const char *ELIXIR_INTERFACE::color_names_table[ELIXIR_INTERFACE::MAX_COLORS] = {
  "black",
  "slateblue",
  "hotpink",
  "lightblue",
  "orange",
  "cyan",
  "blueviolet",
  "yellow",
  "blue",
  "green",
  "blueviolet",
  "red",
  "greenyellow"
}; 

/*
  Standard customization routine.
  Standard Elixir customization routine.  This needs to be defined
  by the calling program to fill the customization menu.
 */
void
ESICustomize(Widget parent_pane)
{
  if (ELIXIR_INTERFACE::client && ELIXIR_INTERFACE::callback)
    (*ELIXIR_INTERFACE::callback)(ELIXIR_INTERFACE::client, parent_pane);
}
#endif
