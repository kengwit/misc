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

#if defined (SILO) && SILO

#include "silo_writer.h"

#undef Min
#define Min(A,B) (((A)<(B))?(A):(B))
#undef Max
#define Max(A,B) (((A)<(B))?(B):(A))

bool
SILO_WRITER::write_db (DBfile *dbfile)
{
  int ndims = 3; // 3D data
#if ! defined(FLT_MAX)
  const double FLT_MAX = 1e300;
#endif  

  // ######################################################
  // Start with the zone list
  int nzones = _nzones;
  int lnodelist = _nodelist.size();
  int *nodelist = new int [lnodelist];
  { int j = 0; DBG( { cerr << "nodelist=[ " << endl; } );
  for (nodelist_vector::iterator i = _nodelist.begin();
       i != _nodelist.end(); i++) {
    nodelist[j] = *i; DBG ( { cerr << *i << ","; } );
    j++;
  } DBG ( { cerr << "]" << endl; } );
  }
  int origin = 1;
  CHECK (_shapecnts.size() == _shapesizes.size(), EXCEPTION_BAD_VALUE,;);
  CHECK (_shapes.size() == _shapesizes.size(), EXCEPTION_BAD_VALUE,;);
  int nshapes = _shapecnts.size();
  int *shapesize = new int [nshapes];
  { int j = 0; DBG( { cerr << "shapesize=[ " << endl; } );
  for (shape_sizes::iterator i = _shapesizes.begin();
       i != _shapesizes.end(); i++) {
    shapesize[j] = *i; DBG ( { cerr << *i << ","; } );
    j++;
  } DBG ( { cerr << "]" << endl; } );
  }
  int *shapecnt = new int [nshapes];
  { int j = 0; DBG( { cerr << "shapecnt=[ " << endl; } );
  for (shape_counts::iterator i = _shapecnts.begin();
       i != _shapecnts.end(); i++) {
    shapecnt[j] = *i; DBG ( { cerr << *i << ","; } );
    j++;
  } DBG ( { cerr << "]" << endl; } );
  }
  DBPutZonelist (dbfile, "zl1", nzones, ndims, nodelist, lnodelist,
                 origin, shapesize, shapecnt, nshapes);
  DBG ( { cerr << "DBPutZonelist done" << endl; } );
  // ######################################################

  // ######################################################
  // Now put the mesh
  char          *coordnames[3];
  coordnames[0] = "xcoords";
  coordnames[1] = "ycoords";
  coordnames[2] = "zcoords";
  
  const int nnodes = _fen_map.size ();
  float         *coords[3];
  float         *x, *y, *z;
  x = new float [nnodes];
  y = new float [nnodes];
  z = new float [nnodes];
  coords[0] = x;
  coords[1] = y;
  coords[2] = z;
  for (fen_map::iterator i = _fen_map.begin(); i != _fen_map.end(); i++) {
    SILO_FEN_DATA_PRIVATE fd = i->second;
    sizet n = fd.private_id - 1; // zero-adjusted
    x[n] = fd.fen_data.x;
    y[n] = fd.fen_data.y;
    z[n] = fd.fen_data.z;
  }
  
  char meshname[] = "famuls_mesh";

  DBoptlist *optlist = DBMakeOptlist (3);
  DBAddOption(optlist, DBOPT_XLABEL, static_cast<void *> ("X Axis"));
  DBAddOption(optlist, DBOPT_YLABEL, static_cast<void *> ("Y Axis"));
  DBAddOption(optlist, DBOPT_ZLABEL, static_cast<void *> ("Z Axis"));
  
  DBPutUcdmesh (dbfile, meshname, ndims, coordnames, coords,
                nnodes, nzones, "zl1", "fl1", DB_FLOAT, optlist);
  DBG ( { cerr << "DBPutUcdmesh done" << endl; } );
  // ######################################################
  
  // ######################################################
  // Material
  char matname[] = "mat1";
  int *matlist = new int [nzones];
  int nmats = 1;
  int matnos[1] = {1};
  int mixlen = 0;
  { for (int j = 0; j < nzones; j++) matlist[j] = matnos[0]; }
  DBPutMaterial (dbfile, matname, meshname, nmats, matnos,
                 matlist, &nzones, 1, NULL, NULL, NULL,
                 NULL, mixlen, DB_FLOAT, optlist);
  DBG ( { cerr << "DBPutMaterial done" << endl; } );
  // ######################################################

  // ######################################################
  // Facelist
  DBfacelist *fl = DBCalcExternalFacelist (nodelist, nnodes, origin,
                                           shapesize, shapecnt, nshapes,
                                           matlist, 0);
  DBG ( { cerr << "DBCalcExternalFacelist done" << endl; } );
  int nfaces = fl->nfaces; DBG ( { cerr << "nfaces " << fl->nfaces << endl; }  );
  int nfshapes = fl->nshapes;
  int *fshapecnt = new int [nfshapes];
  int *fshapesize = new int [nfshapes];
  for (int i = 0; i < nfshapes; i++) {
    fshapecnt[i]  = fl->shapecnt[i];
    fshapesize[i] = fl->shapesize[i];
  }
  int lfacelist = fl->lnodelist;
  int *facelist = new int [lfacelist];
  for (int i = 0; i < lfacelist; i++)
    facelist[i] = fl->nodelist[i];
  int *zoneno = new int [nfaces];
  for (int i = 0; i < nfaces; i++)
    zoneno[i] = fl->zoneno[i];

  DBPutFacelist (dbfile, "fl1", nfaces, ndims, facelist, lfacelist, origin,
                 zoneno, fshapesize, fshapecnt, nfshapes,
                 NULL, NULL, 0); // types, typelist, ntypes
  DBG ( { cerr << "DBPutFacelist done" << endl; } );
  
  DBFreeFacelist(fl);
  // ######################################################

  // ######################################################
  // Put the variables
  for (sizet ivar = 0; ivar < _num_vars; ivar++) {
    char varname[512];
    sprintf (varname, "%s", _var_names[ivar].c_str());
    float *v = new float [nnodes];
    { for (int i = 0; i < nnodes; i++) v[i] = 0; } // zero out
    // and now set to the actual values
    for (fen_map::iterator i = _fen_map.begin(); i != _fen_map.end(); i++) {
      SILO_FEN_DATA_PRIVATE fd = i->second;
      sizet n = fd.private_id - 1; // zero-adjusted
      v[n] = fd.fen_data.varvals[ivar];
    }
    char *varnames[1];
    varnames[0] = varname;
    float *vars[1];
    vars[0] = v;
    DBPutUcdvar (dbfile, varname, meshname, 1, varnames, vars,
                 nnodes, NULL, 0, DB_FLOAT, DB_NODECENT, optlist);
    DBG ( { cerr << "DBPutUcdvar done" << endl; } );
    delete [] v;
  }
  // ######################################################
  
  DBFreeOptlist(optlist);
  
  /*
   * Free the temporary storage.
   */
  delete []  x;
  delete []  y;
  delete []  z;
  delete []  matlist;
  delete [] shapecnt;
  delete [] shapesize;
  delete [] nodelist;
  delete [] fshapecnt;
  delete [] fshapesize;

  return true;
}

static void
report_silo_error (char *msg) {
  cerr << "SILO error: " << msg << endl;
}

void
SILO_WRITER::setup (sizet num_vars, std::vector <std::string> var_names) {
  // The nodal variable arrays are statically allocated:
  // if more is needed, the capacity needs to be increased in the
  // class definition.  (Yes, I know it sucks...)
  CHECK (num_vars < MAX_VARS, EXCEPTION_BAD_VALUE,;);
  // Initialize data structures
  _n_pushed_nodes_per_shape = 0;
  _shapecnts.clear();
  _shapesizes.clear();
  _shapes.clear ();
  _fen_map.clear();
  _nodelist.clear();
  _nzones = 0;
  _curr_shape = ZONE_SHAPE_NONE;
  _num_vars = num_vars;
  // Must get the names for all variables
  CHECK (var_names.size() == _num_vars, EXCEPTION_BAD_VALUE,;);
  _var_names.resize (var_names.size());
  copy (var_names.begin(), var_names.end(), _var_names.begin());
}  

bool
SILO_WRITER::write_silo (std::string filename) {
  // Make sure we're done with the last shape
  if (_curr_shape != ZONE_SHAPE_NONE) end_shape();
  // Set error handling
  DBShowErrors (DB_ALL, report_silo_error);
  // Open the file
  DBfile * dbfile = DBCreate ((char *) filename.c_str(), // file name
                              DB_CLOBBER, // creation mode
                              DB_LOCAL, // destination file format
                              "FAMULS data", // file info
                              DB_PDB // file type
                              );
  // Write it
  write_db (dbfile);
  // Get out
  DBClose(dbfile);
  return true;
}

#endif // Was SILO enabled?

