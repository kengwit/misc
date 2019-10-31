#include "watchpoint.h"
#include "timeP.h"
#include <fstream>
#include <ctime>

WATCH_POINT::WATCH_POINT (list<string> var, string prefix){
  _vbuf.clear();
  for (list<string>::iterator i=var.begin(); i!= var.end() ; i++) {
    class WATCH_POINT::var_buf vbuf;
    vbuf.var = (*i);
    vbuf.first_time = true;
    vbuf.tvp.clear();
    _vbuf.push_back(vbuf);
    _var.push_back(*i);
  }
  _prefix         = prefix;
  _buf_size       = 20;
}

bool WATCH_POINT::set_time_value (string var,double val, double time) {
  for (vector<var_buf>::iterator i = _vbuf.begin(); i != _vbuf.end() ; i++) {
    if ( (*i).var == var) {
      class WATCH_POINT::time_value_pair p;
      p.time  =  time;
      p.value =  val; 
      (*i).tvp.push_back(p);
      if (buffer_size() >= _buf_size) {
        flush (); 
      }
      return 1;
    }    
  }
  return 0;
}

string WATCH_POINT::make_file_name(string var) {
 string filename = _prefix +"_"+ var  + ".xmgr";
 return  filename;
} 

void WATCH_POINT::flush() {
  write_buffer_to_file();
  for (vector<var_buf>::iterator i = _vbuf.begin(); i != _vbuf.end() ; i++) {
    (*i).tvp.clear();
  }
}

void WATCH_POINT::write_buffer_to_file() {
  for(vector<var_buf>::iterator i1 = _vbuf.begin(); i1!= _vbuf.end(); i1++) {
    string var = (*i1).var;
    string file_name = make_file_name(var);
    if((*i1).first_time) {
      write_header (file_name, (*i1).var);
      (*i1).first_time = false;
    } 
    ofstream out; out.open(file_name.c_str(),ios::app);
    for (vector<class WATCH_POINT::time_value_pair>::iterator i=(*i1).tvp.begin(); i!= (*i1).tvp.end(); i++) {
      out<<(*i).time<<" "<<(*i).value<<endl;
    } 
    out.close();
  } 
}

void
WATCH_POINT::write_header (string file, string var)         
{
  int fd;
  ofstream out;

#define POL(aLine) out<<aLine
  out.open(file.c_str(),ios::trunc);
  if(!out.is_open()) {
    string s = "History file can not be emptied\n";
    EXCEPTION_CANNOT_TRUNC e(s);
    throw e;
  }
  out.close();
  out.open(file.c_str(),ios::app);
  if (!out.is_open()){
    string s1 = "Can not append history file\n";
    EXCEPTION_CANNOT_OPEN e1(s1);
    throw e1;
  } 
  
  POL("@    title \"FAMULS (C) 1994-2003 P. Krysl\"\n");
  POL("@g0  autoscale type SPEC\n");
  POL("@    title font 2\n");
  POL("@    title size 1.12\n");
  string date (ckit_get_date_and_time());
  string buffer = "@    subtitle \"" + date + "; " + file + "\"\n";
  POL(buffer);
  POL("@    subtitle font 2\n");
  POL("@    subtitle size 0.8\n");
  POL("@    xaxis  label \"Time [s]\"\n");
  POL("@    xaxis  ticklabel prec 6\n");
  POL("@    xaxis  ticklabel format general\n");
  POL("@    xaxis  ticklabel font 2\n");
  POL("@    xaxis  ticklabel char size 0.67\n");
  POL("@    xaxis  ticklabel start 0.000000\n");
  POL("@    xaxis  tick major grid on\n");
  POL("@    yaxis  label \"" + var +"\"\n");
  POL("@    yaxis  ticklabel prec 6\n");
  POL("@    yaxis  ticklabel format general\n");
  POL("@    yaxis  ticklabel font 2\n");
  POL("@    yaxis  ticklabel char size 0.67\n");
  POL("@    yaxis  ticklabel start 0.000000\n");
  POL("@    yaxis  tick major grid on\n");
  out.close();
}




sizet WATCH_POINT::buffer_size() {
    for (vector<var_buf>::iterator i = _vbuf.begin(); i != _vbuf.end() ; i++) {
      return (*i).tvp.size();  /* all the same size */
    } 
    return 0;
  }
