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
#ifndef CONN_H
# define CONN_H

#include <vector>
#include "conn_base.h"

template <CONN_MANIFOLD_DIM MANIFOLD_DIM, int NFENS, int NREFFENS>
class CONN : public CONN_BASE {

 public: // object functions ////////////////////////////////////////

  /**
   */
  CONN () : CONN_BASE () {}
  ~CONN () {}

  /**
     What is the manifold dimension of the conn?
     (Point=0, line=1, surface=2, or solid=3.) 
   */
  CONN_MANIFOLD_DIM manifold_dim () const { return MANIFOLD_DIM; }

  /**
     Return the number of nodes that this connectivity binds.
   */
  sizet nfens () const { return NFENS; }
  /**
     Return the i-th node.
  */
  FEN *fen (sizet i) const { return _fens(i); }
  /**
     Return the i-th node for modification.
  */
  FEN *& fen (sizet i) { return _fens(i); }
  /**
     Return the local index j: 0 <= j < nfens, of the FEN
     given on input.  A negative number is returned if
     the FEN is not included in the connectivity.
   */
  int local_index (FEN *fen) 
    {
      for (sizet j = 0; j < NFENS; j++) {
        if (fen == _fens(j)) return j;
      }
      return -1;
    }
  /**
     Return the number of refinement nodes that this connectivity uses
     in one refinement step.  See the comments for each template instantiation
     below.
   */
  sizet nreffens () const { return NREFFENS; }
  /**
     Get connectivity code.  This should be a unique number.
     A possibility is to add another template, and set the code in that
     way.  However, the onus would be on the programmer who would have
     to ensure uniqueness.  Current implementation computes the code
     from the three template arguments (admittedly very simplistic).
  */
  CONN_CODE conn_code () { return (((int)MANIFOLD_DIM*1000000) + ((int) NREFFENS * 1000) + NFENS); }

  /**
     Get a refinement node.
     The refinement nodes are oriented for the connectivity on
     input.  Therefore, first try to see if there is any topological
     transformation that may be applied to the connectivity on input
     that would lead to a successful match with self.  If that is
     the case, the same topological operations need to be applied
     to the index into vector of refinement nodes indx, and
     ref_fens[transf(indx)] is returned; otherwise, that is if match
     is not possible, null is returned.
  */
  FEN * get_ref_fen (CONN_BASE *conn, std::vector <FEN *> ref_fens, sizet indx);

  /**
     Clone connectivity.
  */
  CONN_BASE *clone () {
      CONN<MANIFOLD_DIM, NFENS, NREFFENS> *newconn = new CONN<MANIFOLD_DIM, NFENS, NREFFENS> ();
      for (sizet j = 0; j < nfens (); j++)
        newconn->fen(j) = fen(j);
      return newconn;
    }

 private: // object data ////////////////////////////////////////////

  PTR_VECTOR<FEN*, NFENS> _fens;
  
};

/**
   Connectivity of an isolated node.
   The refinement is just another isolated node.
*/
typedef CONN<CONN_0_MANIFOLD, 1, 1> CONN_POINT_1;

/**
   Connectivity of a 2-node line.  
   Numbering of nodes: 
   0-------------------1
   The refinement consists of one interior node, and two vertex nodes.
   X---------0---------X
*/
typedef CONN<CONN_1_MANIFOLD, 2, 1> CONN_LINE_2;

/**
   Connectivity of a 3-node line.  
   Numbering of nodes: 
   0---------2---------1
   The refinement consists of two interior nodes, and three vertex nodes.
   X----0----X----1----X
*/
typedef CONN<CONN_1_MANIFOLD, 3, 2> CONN_LINE_3;

/**
   Connectivity of a 3-node triangle.
   2
   |\
   | \
   |  \
   |   \
   |    \
   0-----1
   The refinement consists of three edge nodes, no interior nodes.
   2
   |\
   | \
   X  X
   |   \
   |    \
   0--X--1
   There are no interior refinement nodes.
*/
typedef CONN<CONN_2_MANIFOLD, 3, 0> CONN_SURF_3;

/**
   Connectivity of a 6-node triangle.
   2
   |\
   | \
   |  \
   5   4
   |    \
   |     \
   |      \
   0---3---1
   There are three interior refinement nodes.
   X
   |\
   | \
   |  \
   X 2 X
   |    \
   | 0 1 \
   |      \
   X---X---X
*/
typedef CONN<CONN_2_MANIFOLD, 6, 3> CONN_SURF_6;

/**
   Connectivity of a 4-node quadrilateral.
   3-------2
   |       |
   |       |
   |       |
   0-------1
   Uses just one interior node for refinement.
   X---X---X
   |       |
   X   0   X
   |       |
   X---X---X
*/
typedef CONN<CONN_2_MANIFOLD, 4, 1> CONN_SURF_4;

/**
   Connectivity of a hexahedron with 8 nodes. 
   Uses just one interior node for refinement.
*/
typedef CONN<CONN_3_MANIFOLD, 8, 1> CONN_SOLID_8;

/**
   Connectivity of a tetrahedron with 4 nodes.
   Does not use any interior node for refinement.
*/
typedef CONN<CONN_3_MANIFOLD, 4, 0> CONN_SOLID_4;

/**
   Connectivity of a tetrahedron with 10 nodes.
   Uses one interior node for refinement.
*/
typedef CONN<CONN_3_MANIFOLD, 10, 1> CONN_SOLID_10;

/**
   Connectivity of a hexahedron with 20 nodes. 
   Uses 7 interior nodes for refinement.
*/
typedef CONN<CONN_3_MANIFOLD, 20, 7> CONN_SOLID_20;

template <>
FEN * CONN_SOLID_10::get_ref_fen (CONN_BASE *conn, vector <FEN *> ref_fens, sizet indx);

template <>
FEN * CONN_SURF_3::get_ref_fen (CONN_BASE *conn, vector <FEN *> ref_fens, sizet indx);

template <>
FEN * CONN_POINT_1::get_ref_fen (CONN_BASE *conn, vector <FEN *> ref_fens, sizet indx);


template <>
FEN * CONN_LINE_2::get_ref_fen (CONN_BASE *conn, vector <FEN *> ref_fens, sizet indx);

template <>
FEN * CONN_LINE_3::get_ref_fen (CONN_BASE *conn, vector <FEN *> ref_fens, sizet indx);

template <>
FEN * CONN_SURF_4::get_ref_fen (CONN_BASE *conn, vector <FEN *> ref_fens, sizet indx);

template <>
FEN * CONN_SURF_6::get_ref_fen (CONN_BASE *conn, vector <FEN *> ref_fens, sizet indx);

template <>
FEN * CONN_SOLID_8::get_ref_fen (CONN_BASE *conn, vector <FEN *> ref_fens, sizet indx);

template <>
FEN * CONN_SOLID_4::get_ref_fen (CONN_BASE *conn, vector <FEN *> ref_fens, sizet indx);

template <>
FEN * CONN_SOLID_20::get_ref_fen (CONN_BASE *conn, vector <FEN *> ref_fens, sizet indx);


#endif
