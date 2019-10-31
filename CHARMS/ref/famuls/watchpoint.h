#ifndef WATCH_POINT_H
# define WATCH_POINT_H

#include "gcell.h"
#include "famexception.h"

class WATCH_POINT {
  
 public: 
  
  WATCH_POINT (list<string> var, string prefix);

  ~WATCH_POINT () {
    _vbuf.clear();
    _var.clear();
  }
  
  list<string> var() {return _var;}
  void set_buf_size (sizet buf_size) { _buf_size = buf_size; }
  bool set_time_value (string var,double val, double time);

  void flush();

 private:
  
  list<string>   _var;
  string         _prefix;
  sizet          _buf_size;
  
 public:
  
  class time_value_pair {
  public:
    double time;
    double value;
  } ;
  class var_buf {
  public:
    bool                      first_time;
    string                    var;
    vector <time_value_pair> tvp;
  } ;
  vector <var_buf> _vbuf; 
  
 private:
  
  string make_file_name(string var);
  void write_buffer_to_file();
  void write_header (string file, string var);
  sizet buffer_size();
  
  DECLARE_EXCEPTION_1_ARG(EXCEPTION_CANNOT_TRUNC, std::string, << arg1);
  DECLARE_EXCEPTION_1_ARG(EXCEPTION_CANNOT_OPEN, std::string, << arg1);
  
};

#endif
