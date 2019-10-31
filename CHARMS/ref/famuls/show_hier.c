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
#include "listP.h"
#define VECMAC_VECTOR_TYPE 1
#include "vecmacP.h" 
#define WANT_VECTP_ALIASES
#include "vectP.h" /* to be able to use vector library routines. */
#include "randP.h" /* to be able to use rand library routines. */
#include "newP.h"              /* Ckit library */
#include "ioP.h"              /* Ckit library */
#include "listP.h"              /* Ckit library */
#include "dictP.h"              /* Ckit library */
#include "getoptP.h"              /* Ckit library */
#include "tokensP.h"              /* Ckit library */
#include "timeP.h"              /* Ckit library */
#include "blpoolP.h"              /* Ckit library */

#  include "Esimple.h"

static int report_opts = 0;
static int nlevels = 8;

static int
usage (int argc, char **argv)
{
  char *option, *usage;

  fprintf (stderr, "usage: %s -i file + options below\n", argv[0]);
  
  if (ckit_opt_first_recognized (&option, &usage)) {
    fprintf (stderr, "%s %s\n",
             (option != NULL ? option : "?"),
             (usage  != NULL ? usage  : "?"));
  }
  while (ckit_opt_next_recognized (&option, &usage)) {
    fprintf (stderr, "%s %s\n",
             (option != NULL ? option : "?"),
             (usage  != NULL ? usage  : "?"));
  }
  exit (0);
  return argc;
}

static void 
process_options (int argc, char **argv)
{
  int help = FALSE;

#define DOSTRINGZ(s) #s
#define STRINGZ(s) DOSTRINGZ(s)
#define SET_BOOL(name)                                          \
  {                                                             \
    if (report_opts) { fprintf (stderr, "Set: %s\n", #name); }  \
    name = TRUE;                                                \
  }
#define SET_WITH_ARG(name, format)                                        \
  {                                                                       \
    if (report_opts) { fprintf (stderr, "Set: %s = %" #format "\n",       \
                                STRINGZ(name), name); }                   \
  }
  
  GET_OPT_BOOL (argc, argv, "-report_opts", SET_BOOL(report_opts),
                "prints options as they are being set");
  
  GET_OPT_BOOL (argc, argv, "-help", SET_BOOL(help),
                "prints usage for each option recognized");

  nlevels = 8;
  GET_OPT_ARG (argc, argv, "-nlevels", int, nlevels,
               { SET_WITH_ARG (nlevels, d); },
               "number of levels");
  
  if (help)
    usage (argc, argv);
}

static void
toggle_layer_cb (Widget w, XtPointer ptr, XtPointer call_data)
{
  Boolean s;
  Arg al[1];
  int ac;
  char buf[132];
  int layer = (int)ptr;
  
  ac = 0;
  XtSetArg(al[ac], XtNstate, &s); ac++;
  XtGetValues(w, al, ac);
  if (s)
    sprintf(buf, "LAYER ON %d", layer);
  else
    sprintf(buf, "LAYER OFF %d", layer);
  ESIHandleCmd(buf);
}

static void
animate_cb (Widget w, XtPointer ptr, XtPointer call_data)
{
  static int snapshot = 0;
  int i;
  const int maxi = 20;
  char buf[128];
  for (i = 0; i < maxi; i++) {
    ESIHandleCmd ("view rotnor z 0.01");
    ESIHandleCmd ("view redraw");
    sprintf (buf, "image ss%05d.jpg", snapshot++);
    sleep (2);
    ESIHandleCmd (buf);
  }
}

static void
zoom_in_cb (Widget w, XtPointer ptr, XtPointer call_data)
{
  ESIHandleCmd ("view zoom 0.95");
}

static void
zoom_out_cb (Widget w, XtPointer ptr, XtPointer call_data)
{
  ESIHandleCmd ("view zoom 1.05");
}

void
ESICustomize (Widget parent_pane)
{
  int j, argn;
  Arg args[13];
  char buf[32];
  Widget level_visibility_palette;

  argn = 0;
  ESIAddPalette("level_visibility", "Level visibility ...", parent_pane,
                args, argn, &level_visibility_palette);
  
  for (j = 1; j <= nlevels; j++) {
    sprintf (buf, "Level %3d", j);
    argn = 0;
    XtSetArg (args[argn], XtNstate, TRUE); argn++;
    ESIAddButton("level_visibility_toggle", buf, 
                 toggleWidgetClass, level_visibility_palette,
                 args, argn, toggle_layer_cb, (XtPointer)j);
  }
  
  argn = 0;
  ESIAddButton("animate",
                   "Animate",
                   commandWidgetClass, parent_pane,
                   NULL, 0, animate_cb, (XtPointer)NULL);
  argn = 0;
  ESIAddButton("zoom_out",
                   "Zoom out",
                   commandWidgetClass, parent_pane,
                   NULL, 0, zoom_out_cb, (XtPointer)NULL);
  argn = 0;
  ESIAddButton("zoom_in",
                   "Zoom in",
                   commandWidgetClass, parent_pane,
                   NULL, 0, zoom_in_cb, (XtPointer)NULL);
}

static int
has_suffix (char *f, char *s /* suffix */)
{
  int lf, ls, i;

  if (!f || !s)
    return FALSE;
  
  lf = strlen((f)); ls = min(strlen((s)), lf);
  for (i = 0; i < ls; i++)
    if (toupper((f)[lf-i]) != toupper((s)[ls-i]))
      return FALSE;
  return TRUE;
}

static void
suck_in_file (char *f)
{
  char command[1024];
  
  if (has_suffix (f, ".egf")) {
    sprintf (command, "LOAD \"%s\"", f);
    ESIHandleCmd (command);
  }  else
    fprintf (stderr, "Did not recognize format of \"%s\"\n", f);
}

static void
process_args (int argc, char **argv)
{
  int i;
  i = 1;
  while (i < argc) {
    if (ckit_io_file_exists (argv[i]))
      suck_in_file (argv[i]);
    i++;
  }
}

int
main (int argc, char **argv)
{
  unsigned long m = ESI_GRAPHIC_EDITOR_MASK;

  process_options (argc, argv);

  ESIBuildInterface (m, argc, argv);
  ESIPopup ();
  process_args (argc, argv);
  ESIHandleCmd ("EDGECOLOR white");
  ESIHandleCmd ("EDGEFLAG on");
  ESIHandleCmd ("COLOR yellow");
  ESIHandleCmd ("MSIZE 10");
  ESIHandleCmd ("MTYPE 1");
  ESIHandleCmd ("FONT \"9x15\"");
  ESIHandleCmd ("SELCRIT inside");
  ESIHandleCmd ("newframe width 480 height 480 x 400 y 10");
  ESIHandleCmd ("view bground grey50");
  ESIHandleCmd ("view axes off");
  ESIHandleCmd ("fit");
  ESIPopupAndRun ();
  
  return 0;
}
