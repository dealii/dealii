#//---------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------


#include <base/vector_slice.h>
#include <base/logstream.h>

#include <vector>
#include <fstream>
#include <iostream>

void f(const std::vector<int>& v)
{
  const VectorSlice<const std::vector<int> >
    s = make_slice(v,2,3);
  
  for (unsigned int i=0;i<s.size();++i)
    std::cerr << '\t' << s[i];
  std::cerr << std::endl;
}


int main()
{
  std::ofstream logfile("slice_vector.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

                                   // we do something rather weird in
                                   // this file: bind the buffer of
                                   // cerr to the log file, so that
                                   // the output of the Assert* macros
                                   // that are written to std::cerr
                                   // end up in the logfiles.
                                   //
                                   // so make sure we store a pointer
                                   // to the old buffer, switch to the
                                   // buffer of the logfile, and at
                                   // the end of main() switch
                                   // back. note that if we don't
                                   // switch back, we get a segfault
                                   // later on in the destruction of
                                   // std::cerr, since it tries to do
                                   // something with its buffer, but
                                   // that is already gone by then
#if __GNUC__ != 2
  std::basic_streambuf<char> *old_cerr_buf = std::cerr.rdbuf();
#else
  streambuf *old_cerr_buf = std::cerr.rdbuf();
#endif
  std::cerr.rdbuf(logfile.rdbuf());
  
  std::vector<int> v(7);

  for (unsigned int i=0;i<v.size();++i)
    v[i] = i;
  
  VectorSlice<std::vector<int> > s(v, 3, 4);

  for (unsigned int i=0;i<s.size();++i)
    s[i] = i;

  for (unsigned int i=0;i<v.size();++i)
    std::cerr << '\t' << v[i];
  std::cerr << std::endl;

  f(v);

  const VectorSlice<const std::vector<int> >
    s2 = make_slice(v, 3, 5);
  int n = s[4];
  n += 3;

  std::cerr.rdbuf(old_cerr_buf);
}
