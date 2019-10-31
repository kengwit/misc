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
#ifndef UTIL_H
# define UTIL_H

#define DO_STRINGIZE_(S) #S
#define STRINGIZE_(S) DO_STRINGIZE_(S)
#define EXP(S) STRINGIZE_(S)

/* Ordinary loop */
#define FOR(_Counter_, _Up_) for (_Counter_ = 0; _Counter_ < _Up_; _Counter_++)

#define WARN(msg) { cerr << "**** WARNING (" << __FILE__ << " l. " << __LINE__ << ") ****: " << msg << endl; }
#define WARN2(msg,arg2) { cerr << "**** WARNING (" << __FILE__ << " l. " << __LINE__ << ") ****: " << msg << " " << arg2 << endl; }

#if !defined(max)
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#if !defined(min)
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#endif
