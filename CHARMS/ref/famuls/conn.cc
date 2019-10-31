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
#include "conn.h"

using namespace std;

#undef SWAP
#define SWAP(n0, n1) { FEN *afen = fen(n0); fen(n0) = fen(n1); fen(n1) = afen; }

template <>
FEN * CONN_POINT_1::get_ref_fen (CONN_BASE *conn, vector <FEN *> ref_fens, sizet indx)
{
  CHECK (ref_fens.size () == 1, EXCEPTION_BAD_VALUE,;);
  CHECK (conn->manifold_dim () == manifold_dim (), EXCEPTION_BAD_VALUE,;);
  CHECK (conn->nfens () == nfens (), EXCEPTION_BAD_VALUE,;);
  CHECK (indx == 0, EXCEPTION_BAD_VALUE ,;);
  if (conn->fen(0) == fen(0)) return ref_fens[0];
  else                        return 0;
}


template <>
FEN * CONN_LINE_2::get_ref_fen (CONN_BASE *conn, vector <FEN *> ref_fens, sizet indx)
{
  CHECK (ref_fens.size () == 1, EXCEPTION_BAD_VALUE,;);
  CHECK (conn->manifold_dim () == manifold_dim (), EXCEPTION_BAD_VALUE ,;);
  CHECK (conn->nfens () == nfens (), EXCEPTION_BAD_VALUE ,;);
  CHECK (indx == 0, EXCEPTION_BAD_VALUE ,;);
  if        (conn->fen(0) == fen(0) && conn->fen(1) == fen(1)) {
    return ref_fens[0];
  } else if (conn->fen(0) == fen(1) && conn->fen(1) == fen(0)) {
    return ref_fens[0];
  } else {
    return 0;
  }
}

template <>
FEN * CONN_LINE_3::get_ref_fen (CONN_BASE *conn, vector <FEN *> ref_fens, sizet indx)
{
  CHECK (ref_fens.size () == 2, EXCEPTION_BAD_VALUE ,;);
  CHECK (conn->manifold_dim () == manifold_dim (), EXCEPTION_BAD_VALUE ,;);
  CHECK (conn->nfens () == nfens (), EXCEPTION_BAD_VALUE ,;);
  CHECK ((indx == 0) || (indx == 1), EXCEPTION_BAD_VALUE ,;);
  if        (conn->fen(0) == fen(0) && conn->fen(1) == fen(1) && conn->fen(2) == fen(2)) {
    return ref_fens[indx];
  } else if (conn->fen(0) == fen(1) && conn->fen(1) == fen(0) && conn->fen(2) == fen(2)) {
    return ref_fens[(indx+1)%2]; // reverse
  } 
  return 0;
}

template <>
FEN * CONN_SURF_3::get_ref_fen (CONN_BASE *conn, vector <FEN *> ref_fens, sizet indx)
{
  CHECK (ref_fens.size () == 1, EXCEPTION_BAD_VALUE ,;);
  CHECK (conn->manifold_dim () == manifold_dim (), EXCEPTION_BAD_VALUE ,;);
  CHECK (conn->nfens () == nfens (), EXCEPTION_BAD_VALUE ,;);
  return 0; // The triangle has no internal refinement nodes
}

template <>
FEN * CONN_SURF_4::get_ref_fen (CONN_BASE *conn, vector <FEN *> ref_fens, sizet indx)
{
  CHECK (ref_fens.size () == 1, EXCEPTION_BAD_VALUE ,;);
  CHECK (conn->manifold_dim () == manifold_dim (), EXCEPTION_BAD_VALUE ,;);
  CHECK (conn->nfens () == nfens (), EXCEPTION_BAD_VALUE ,;);
  CHECK (indx == 0, EXCEPTION_BAD_VALUE ,;);
  // Check by trying the four configurations obtained by rotation
#define CHK4(n0,n1,n2,n3) ((conn->fen(0) == fen(n0)) && \
                           (conn->fen(1) == fen(n1)) && \
                           (conn->fen(2) == fen(n2)) && \
                           (conn->fen(3) == fen(n3)))
  if        (CHK4(0,1,2,3)) {
    return ref_fens[0];
  } else if (CHK4(1,2,3,0)) {
    return ref_fens[0];
  } else if (CHK4(2,3,0,1)) {
    return ref_fens[0];
  } else if (CHK4(3,0,1,2)) {
    return ref_fens[0];
  } else if (CHK4(3,2,1,0)) { // and then try the mirror configurations
    return ref_fens[0];
  } else if (CHK4(2,1,0,3)) {
    return ref_fens[0];
  } else if (CHK4(1,0,3,2)) {
    return ref_fens[0];
  } else if (CHK4(0,3,2,1)) {
    return ref_fens[0];
  }
  return 0; // no match
}

template <>
FEN * CONN_SURF_6::get_ref_fen (CONN_BASE *conn, vector <FEN *> ref_fens, sizet indx)
{
  CHECK (ref_fens.size () == 3, EXCEPTION_BAD_VALUE ,;);
  CHECK (conn->manifold_dim () == manifold_dim (), EXCEPTION_BAD_VALUE ,;);
  CHECK (conn->nfens () == nfens (), EXCEPTION_BAD_VALUE ,;);
  CHECK ((indx == 0) || (indx == 1) || (indx == 2), EXCEPTION_BAD_VALUE ,;);
  // Check by trying the four configurations obtained by rotation
#undef CHK6
#define CHK6(n0,n1,n2,n3,n4,n5) ((conn->fen(0) == fen(n0)) &&   \
                                 (conn->fen(1) == fen(n1)) &&   \
                                 (conn->fen(2) == fen(n2)) &&   \
                                 (conn->fen(3) == fen(n3)) &&   \
                                 (conn->fen(4) == fen(n4)) &&   \
                                 (conn->fen(5) == fen(n5))      \
                                 )
  if        (CHK6(0,1,2,3,4,5)) {
    return ref_fens[indx]; // 0 1 2
  } else if (CHK6(1,2,0,4,5,3)) {
    return ref_fens[(indx+1)%3]; // 1 2 0
  } else if (CHK6(2,0,1,5,3,4)) {
    return ref_fens[(indx+2)%3]; // 2 0 1
  }
  // Now the mirror reflections
  if        (CHK6(0,2,1,5,4,3)) {
    const int i[3] = {0,2,1};
    return ref_fens[i[indx]];
  } else if (CHK6(2,1,0,4,3,5)) {
    const int i[3] = {2,1,0};
    return ref_fens[i[indx]]; 
  } else if (CHK6(1,0,2,3,5,4)) {
    const int i[3] = {1,0,2};
    return ref_fens[i[indx]]; 
  }
  return 0; // no match
}

template <>
FEN * CONN_SOLID_8::get_ref_fen (CONN_BASE *conn, vector <FEN *> ref_fens, sizet indx)
{
  CHECK (ref_fens.size () == 1, EXCEPTION_BAD_VALUE ,;);
  CHECK (conn->manifold_dim () == manifold_dim (), EXCEPTION_BAD_VALUE ,;);
  CHECK (conn->nfens () == nfens (), EXCEPTION_BAD_VALUE ,;);
  CHECK (indx == 0, EXCEPTION_BAD_VALUE ,;);
  // We could proceed as with the quadrilateral, but another option
  // is just to count the number of present nodes
  sizet num_matched = 0;
  for (sizet j = 0; j < nfens (); j++) {
    const FEN *to_match = fen(j);
    for (sizet k = 0; k < nfens (); k++) {
      if (to_match == conn->fen(k)) {
        num_matched++;
        break;
      }
    }
  }
  if (num_matched == 8) return ref_fens[0];
  else                  return 0;
}


template <>
FEN * CONN_SOLID_4::get_ref_fen (CONN_BASE *conn, vector <FEN *> ref_fens, sizet indx)
{
  CHECK (ref_fens.size () == 1, EXCEPTION_BAD_VALUE ,;);
  CHECK (conn->manifold_dim () == manifold_dim (), EXCEPTION_BAD_VALUE ,;);
  CHECK (conn->nfens () == nfens (), EXCEPTION_BAD_VALUE ,;);
  return 0; // The 4-node tet has no internal refinement nodes
}

template <>
FEN * CONN_SOLID_10::get_ref_fen (CONN_BASE *conn, vector <FEN *> ref_fens, sizet indx)
{
  CHECK (ref_fens.size () == 1, EXCEPTION_BAD_VALUE ,;);
  CHECK (conn->manifold_dim () == manifold_dim (), EXCEPTION_BAD_VALUE ,;);
  CHECK (conn->nfens () == nfens (), EXCEPTION_BAD_VALUE ,;);
  CHECK (indx == 0, EXCEPTION_BAD_VALUE ,;);
  // We could proceed as with the quadrilateral, but another option
  // is just to count the number of present nodes
  sizet num_matched = 0;
  for (sizet j = 0; j < nfens (); j++) {
    const FEN *to_match = fen(j);
    for (sizet k = 0; k < nfens (); k++) {
      if (to_match == conn->fen(k)) {
        num_matched++;
        break;
      }
    }
  }
  if (num_matched == 10) return ref_fens[0];
  else                   return 0;
}


template <>
FEN * CONN_SOLID_20::get_ref_fen (CONN_BASE *conn, vector <FEN *> ref_fens, sizet indx)
{
  CHECK (indx == 0, EXCEPTION_NOT_IMPLEMENTED,;);
  return 0;
}

