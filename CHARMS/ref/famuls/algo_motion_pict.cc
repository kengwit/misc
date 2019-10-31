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
#include "algo_motion_pict.h"
#include "args.h"
#include "field.h"
#include "evalpt.h"

const double ALGO_MOTION_PICT::SHRINK = 0.9;
#if defined (GIFACE) && GIFACE
bool ALGO_MOTION_PICT::initialized = false;

static void 
set_displ_scale(Widget w, XtPointer client_ptr, XtPointer call_data)
{
  if (TypeInParseLine((char *) call_data) > 0 && TypeInGetTokenType(1) == NUMBER) {
    ALGO_MOTION_PICT *a = (ALGO_MOTION_PICT *) EMGetClient (ESIModel());
    a->set_scale(TypeInGetTokenAsDouble(1));
  }
}

static void 
StopComputation(Widget w, XtPointer text_ptr, XtPointer call_data)
{
  ESIEventLoop(1, "Stopped: press ctrl-p or push button Proceed");
} 


static void 
Proceed(Widget w, XtPointer text_ptr, XtPointer call_data)
{
  ESIEventLoopProceed();
}

static void
redisplay_CB(Widget w, XtPointer ptr, XtPointer call_data)
{
}


/*
  Standard customization routine.
  Standard Elixir customization routine.  This needs to be defined
  by the calling program to fill the customization menu.
 */
static void
myESICustomize(void *client, Widget parent_pane
             )
{
  Widget settings, graphic_control, probe, modr;
  Arg args[13];
  int argn;
  Boolean b;

  argn = 0;
  ESIAddButton("proceed", "Proceed", commandWidgetClass,
                   parent_pane, NULL, 0, Proceed, NULL);
  argn = 0;
  ESIAddButton("stop", "Stop", commandWidgetClass,
                   parent_pane, NULL, 0, StopComputation, NULL);
  argn = 0;
  ESIAddPopupDialog("displ_scale", "Displacement scaling",
                    "Input the displacement scale",
                    "1.0", parent_pane, args, argn,
                    set_displ_scale, NULL,
                    ESIDialogValueNumber, NULL);
  ESIAddButton("redisplay",
                   "Redisplay",
                   commandWidgetClass, parent_pane,
                   NULL, 0, redisplay_CB, (XtPointer)NULL);
}
#endif

ALGO_MOTION_PICT::ALGO_MOTION_PICT(string name, MGR *mgr, GMESH *gmesh) : ALGO (name, mgr)
{
#if defined (GIFACE) && GIFACE
  if (!initialized) {
    ELIXIR_INTERFACE::build_interface (this, myESICustomize);
    ESIHandleCmd ("color red");
    ESIHandleCmd ("shrink 1.0");
    ESIHandleCmd ("fill solid");
    ESIHandleCmd ("msize 6");
    ESIHandleCmd ("mtype 1");
    ESIHandleCmd ("newframe width 480 height 480 x 10 y 20");
    EMSetClient (ESIModel (), this);
    BOOLEAN suc;
    fen_color = ColorGetPixelFromString ("red", &suc);
    CHECK (suc, EXCEPTION_NULL_PTR,;);
    gcell_color = ColorGetPixelFromString ("red", &suc);
    CHECK (suc, EXCEPTION_NULL_PTR,;);
    edge_color = ColorGetPixelFromString ("white", &suc);
    CHECK (suc, EXCEPTION_NULL_PTR,;);
    gcell_layer = 0;
    fen_layer = 1;
    _fen_rep = FEN::MARKER;
    _node_size_mult = 0.06;
    _node_size = node_size(gmesh);
    _level_as_fen_layer = true;
    _gogl.clear();
    _gofl.clear();
    ESIPopup ();
    initialized = true;
  }
#endif
}

ALGO_MOTION_PICT::~ALGO_MOTION_PICT () {
#if defined (GIFACE) && GIFACE
  _gogl.clear();
  _gofl.clear();
#endif
}

void ALGO_MOTION_PICT::event_loop (bool stop, string msg) {
#if defined (GIFACE) && GIFACE
  ESIEventLoop(stop, (char *) msg.c_str());
#endif
}


#if defined (GIFACE) && GIFACE
double
ALGO_MOTION_PICT::node_size (GMESH *gmesh)
{
  G3D_box_t mesh_box;
  gmesh->box (&mesh_box);
  double r[3];
  r[0] = G3D_BOX_RANGE (mesh_box, x);
  r[1] = G3D_BOX_RANGE (mesh_box, y);
  r[2] = G3D_BOX_RANGE (mesh_box, z);
  sizet fen_total = gmesh->fen_total();
  double mesh_vol = 1;
  int p = 0;
  for (int i = 0; i < 3; i++) {
    if (r[i] > 0) {
      mesh_vol *= r[i];
      p++;
    }
  }
  double ns = pow (mesh_vol / (fen_total + 1), 1.0/p) * _node_size_mult;
  return ns;
}
#endif // GIFACE

void
ALGO_MOTION_PICT::add_graphics (GCELL *gc, FIELD_VECTOR *geometry)
{
#if defined (GIFACE) && GIFACE
  vector <FIXED_VECTOR<3> > p; p.resize(gc->conn()->nfens());
  for (sizet i = 0; i < gc->conn()->nfens(); i++) {
    p[i] = geometry->evaluate(geometry, gc->conn()->fen(i), gc);
  }
  go_gcell_pair_t pr;
  pr.go = gc->make_elixir_obj (p);
  pr.gcell = gc;
  if (pr.go != 0) {
    EASValsSetLayer (gcell_layer); 
    EASValsSetFillStyle (FILL_SOLID);
    EASValsSetColor (gcell_color);
    EASValsSetShrink (SHRINK);
    EASValsSetEdgeColor (edge_color);
    EASValsSetEdgeFlag (TRUE);
    EGWithMaskChangeAttributes (ALL_ATTRIB_MASK, pr.go);
    EMNoDrawAddGraphicsToModel (ESIModel(), pr.go);
    _gogl.push_back(pr);
  }
#endif // GIFACE
}

void
ALGO_MOTION_PICT::add_graphics (FEN *fen, FIELD_VECTOR *geometry)
{
#if defined (GIFACE) && GIFACE
  FIXED_VECTOR<3> p;
  p = geometry->evaluate(geometry, fen);
  go_fen_pair_t pr;
  pr.go = fen->make_elixir_obj (p, _node_size, _fen_rep);
  pr.fen = fen;
  if (pr.go != 0) {
    EASValsSetLayer (fen_layer); 
    EASValsSetFillStyle (FILL_SOLID);
    EASValsSetColor (fen_color);
    EGWithMaskChangeAttributes (ALL_ATTRIB_MASK, pr.go);
    EMNoDrawAddGraphicsToModel (ESIModel(), pr.go);
    _gofl.push_back(pr);
  }
#endif // GIFACE
}

static int 
redraw_view (NODE ptr, NODE elem)
{
#if defined (GIFACE) && GIFACE
  EVFastRedraw ((EView *)elem);
#endif // GIFACE
  return 1;
}

void
ALGO_MOTION_PICT::update (FIELD_VECTOR *geometry)
{
#if defined (GIFACE) && GIFACE
  for (list <go_gcell_pair_t >::iterator it = _gogl.begin();
       it != _gogl.end(); it++) {
    go_gcell_pair_t pr = *it;
    const sizet maxp = 1000;
    WCRec points[maxp]; // RAW
    CHECK(maxp > pr.gcell->conn()->nfens(), EXCEPTION_ILLEGAL_USE,;);
    vector <FIXED_VECTOR<3> > p; p.resize(pr.gcell->conn()->nfens());
    for (sizet j = 0; j < pr.gcell->conn()->nfens(); j++) {
      FIXED_VECTOR<3> p = geometry->evaluate(geometry, pr.gcell->conn()->fen(j), pr.gcell);
      points[j].x = p(0);
      points[j].y = p(1);
      points[j].z = p(2);
    }
    EGModifyGraphicsGeometry(pr.go, points);
  }
  for (list <go_fen_pair_t >::iterator it = _gofl.begin();
       it != _gofl.end(); it++) {
    go_fen_pair_t pr = *it;
    WCRec points[1]; 
    FIXED_VECTOR<3> p = geometry->evaluate(geometry, pr.fen);
    points[0].x = p(0);
    points[0].y = p(1);
    points[0].z = p(2);
    EGModifyGraphicsGeometry(pr.go, points);
  }
  EMDispatchToDependentViews (ESIModel (), redraw_view, NULL);
  event_loop(false," ");
#endif // GIFACE
}
