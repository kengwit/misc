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
#ifndef EXCEPTION_H
#  define EXCEPTION_H


#include <iostream>
#include <exception>
#include <cstdlib>
#include <string>
#include <strstream>
#include <exception>

#ifndef __GNUC__
#  define __FUNCTION__ "(unknown function)"
#endif

class BASE_EXCEPTION : public std::exception {

  public:
    BASE_EXCEPTION ();
    BASE_EXCEPTION (const char* file, const int line, const char *function,
                    const char* condition, const char *exception_name);
    void set_fields (const char *f,
                             const int   l,
                             const char *func,
                             const char *c,
                             const char *e);
    void print_exception_data (std::ostream &out) const;
    virtual void print_info (std::ostream &out) const;
    ~BASE_EXCEPTION () throw() {}

  protected:
    const char  *_file;
    unsigned int _line;
    const char  *_function;
    const char  *_condition;
    const char  *_exception_name;
};


# define CHECK(cond, EXC, ARGS)                                         \
  {                                                                     \
    if (!(cond)) {                                                      \
      EXC e ARGS;                                                       \
      e.set_fields (__FILE__, __LINE__, __FUNCTION__, #cond, #EXC);     \
      std::cerr << "!<<<<<<<<<<" << endl;                               \
      e.print_exception_data (std::cerr);                               \
      e.print_info (std::cerr);                                         \
      std::cerr << "!>>>>>>>>>" << endl;                                \
      abort ();                                                         \
    }                                                                   \
  }


#define CHECK_THROW(cond, EXC, ARGS)                                    \
  {                                                                     \
    if (!(cond)) {                                                      \
      EXC e ARGS;                                                       \
      e.set_fields (__FILE__, __LINE__, __FUNCTION__, #cond, #EXC);     \
      e.print_exception_data (cerr);                                    \
      throw e;                                                          \
    }                                                                   \
  }


#define DECLARE_EXCEPTION_1_ARG(EXC, TYPE1, OUT)                \
   class EXC : public BASE_EXCEPTION {                          \
     public:                                                    \
      EXC (const TYPE1 A1) : arg1 (A1) {};                      \
      virtual void print_info (std::ostream &out) const {       \
        out OUT << std::endl;                                   \
      };                                                        \
    ~EXC () throw() {}                                          \
     private:                                                   \
      const TYPE1 arg1;                                         \
   }

#include "util.h"

namespace STANDARD_EXCEPTIONS 
{
  class EXCEPTION_NOT_IMPLEMENTED :  public BASE_EXCEPTION {};
  class EXCEPTION_NULL_PTR        :  public BASE_EXCEPTION {};
  class EXCEPTION_BAD_ACCESS      :  public BASE_EXCEPTION {};
  class EXCEPTION_BAD_SYNTAX      :  public BASE_EXCEPTION {};
  class EXCEPTION_BAD_VALUE       :  public BASE_EXCEPTION {};
  class EXCEPTION_ILLEGAL_USE     :  public BASE_EXCEPTION {};

  DECLARE_EXCEPTION_1_ARG (EXCEPTION_FILE_ERROR, std::string, << arg1);
  DECLARE_EXCEPTION_1_ARG (EXCEPTION_IN_LOCAL_LIBS, std::string, << arg1);

}

using namespace STANDARD_EXCEPTIONS;

#endif
