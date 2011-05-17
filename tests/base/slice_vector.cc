#//---------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005, 2006, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------


#include "../tests.h"
#include <deal.II/base/vector_slice.h>
#include <deal.II/base/logstream.h>

#include <vector>
#include <fstream>
#include <iomanip>

void f(const std::vector<int>& v)
{
  const VectorSlice<const std::vector<int> >
    s = make_slice(v,2,3);
  
  for (unsigned int i=0;i<s.size();++i)
    deallog << '\t' << s[i];
  deallog << std::endl;
}


int main()
{
  deal_II_exceptions::disable_abort_on_exception();
  std::ofstream logfile("slice_vector/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  std::vector<int> v(7);

  for (unsigned int i=0;i<v.size();++i)
    v[i] = i;
  
  VectorSlice<std::vector<int> > s(v, 3, 4);

  for (unsigned int i=0;i<s.size();++i)
    s[i] = i;

  for (unsigned int i=0;i<v.size();++i)
    deallog << '\t' << v[i];
  deallog << std::endl;

  f(v);

  make_slice(v, 3, 5);
  int n = s[4];
  n += 3;
}
