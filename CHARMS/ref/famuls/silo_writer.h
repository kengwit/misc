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
#ifndef SILO_WRITER_H
#   define SILO_WRITER_H

// Is SILO enabled?
#if defined (SILO) && SILO

#include "famuls.h"
extern "C" {
#include <silo.h>
}
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <vector.h>
#include <map.h>
#include <set.h>
#include "famexception.h"

#define DBG(code) { code; }
/**
   Class for writing SILO post-processing files.
   Protocol:
   o Constructor.
   o Call setup().
   o Use PROTO_SILO to stream the zone data (each ecell will call
     start_zone() and end_zone() to indicate where the streaming of
     node data begins and ends, and .
   o Call write_silo().
   Also see ALGO_SILO.
*/
class SILO_WRITER {

 public: // class definitions //////////////////////////////////

  static const sizet MAX_VARS = 3;
  
  /**
     Data type to transmit the node data to the SILO_WRITER.
  */
  typedef struct {
    sizet id;         // set the id
    float x, y, z;     // set the coordinates
    float varvals[MAX_VARS]; // set the variable values
  } SILO_FEN_DATA;
  
  typedef enum SILO_UCD_ZONE_SHAPE {
    ZONE_SHAPE_NONE = 0,
    // 1D
    ZONE_SHAPE_BEAM = DB_ZONETYPE_BEAM, 
    // 2D
    ZONE_SHAPE_POLYGON = DB_ZONETYPE_POLYGON,
    ZONE_SHAPE_TRIANGLE = DB_ZONETYPE_TRIANGLE,
    ZONE_SHAPE_QUAD = DB_ZONETYPE_QUAD,
    // 3D
    ZONE_SHAPE_POLYHEDRON = DB_ZONETYPE_POLYHEDRON,
    ZONE_SHAPE_TET = DB_ZONETYPE_TET,
    ZONE_SHAPE_PYRAMID = DB_ZONETYPE_PYRAMID,
    ZONE_SHAPE_PRISM = DB_ZONETYPE_PRISM,
    ZONE_SHAPE_HEX = DB_ZONETYPE_HEX
  } SILO_UCD_ZONE_SHAPE;
  
 public: // object methods //////////////////////////////////////////////
  
  /**
     Writer of SILO post-processing files.
   */
  SILO_WRITER () {
    _n_pushed_nodes_per_shape = 0;
    _shapecnts.clear();
    _shapesizes.clear();
    _shapes.clear ();
    _fen_map.clear();
    _nodelist.clear();
    _nzones = 0;
    _curr_shape = ZONE_SHAPE_NONE;
    _num_vars = 0;
    _var_names.clear();
  }
  
  /**
     Use before streaming any data to the writer.
     Supply the number of variables to associate with each
     node.
   */
  void setup (sizet num_vars, std::vector <std::string> var_names);
  
  /**
     Write the data that the writer accumulated so far.
   */
  bool write_silo (std::string filename);

  /**
     Start a zone (zone=ecell in SILO parlance).
   */
  void start_zone (SILO_UCD_ZONE_SHAPE shape) {
    if ((shape != _curr_shape) || (shape == ZONE_SHAPE_POLYHEDRON)) {
      start_shape (shape);
    }
    sizet n = _shapecnts.size()  - 1;
    _shapecnts[n] = _shapecnts[n]+1; // increment number of zones of this shape
    _n_pushed_nodes_per_shape = 0;
  }
  /**
     End a zone.
   */
  void end_zone () {
    _nzones++;
  }

  /**
     Push an integer value to node list.  This
     may be used to stream ZONE_SHAPE_POLYHEDRON data.
   */
  void push_to_nodelist (int i) {
    _n_pushed_nodes_per_shape++;
    _nodelist.push_back (i); 
  }

  /**
     Push node data to the node list.
   */
  void push_to_nodelist (SILO_FEN_DATA *fen_data) {
    SILO_FEN_DATA_PRIVATE fd;
    fd.fen_data = *fen_data;
    _n_pushed_nodes_per_shape++;
    fen_map::iterator i = _fen_map.find (fen_data->id);
    if (i == _fen_map.end()) {
      const sizet n = _fen_map.size();
      fd.private_id = n+1;
      _fen_map.insert (fen_map::value_type (fd.fen_data.id, fd));
    } else {
      fd = i->second;
    }
    _nodelist.push_back (fd.private_id);
    DBG ( { cerr << "pushing " << fd.fen_data.id << " (private " << fd.private_id << ")" << endl; } );
  }
  
 private: // object data ////////////////////////////////////////////////
  
  typedef struct {
    SILO_FEN_DATA fen_data;
    sizet private_id;
    float varvals[MAX_VARS];
  } SILO_FEN_DATA_PRIVATE;

  typedef std::vector <int> shape_counts;
  typedef std::vector <int> shape_sizes;
  typedef std::vector <int> nodelist_vector;
  typedef std::vector <SILO_UCD_ZONE_SHAPE> shape_types;
  typedef std::map <sizet, SILO_FEN_DATA_PRIVATE> fen_map;

  int                        _n_pushed_nodes_per_shape;
  shape_counts               _shapecnts;
  shape_sizes                _shapesizes;
  shape_types                _shapes;
  fen_map                    _fen_map;
  nodelist_vector            _nodelist;
  sizet                     _nzones;
  SILO_UCD_ZONE_SHAPE        _curr_shape;
  sizet                     _num_vars;
  std::vector <std::string>  _var_names;

 private: // object private methods ////////////////////////////////

  bool write_db (DBfile *dbfile);

  /**
     When the shape of the zone that is going to be started is
     different from the shape of the previous zone, start a new
     shape.
   */
  void start_shape (SILO_UCD_ZONE_SHAPE shape) {
    if (_curr_shape != ZONE_SHAPE_NONE) end_shape();
    _shapes.push_back (shape);
    _shapecnts.push_back(0);
    _curr_shape = shape;
  }

  /**
   */
  void end_shape () {
    _shapesizes.push_back (_n_pushed_nodes_per_shape);
    _curr_shape = ZONE_SHAPE_NONE;
  }

};

// Is SILO enabled?
#endif

#endif

