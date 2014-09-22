// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include "../tests.h"
#include <deal.II/base/vector_slice.h>
#include <deal.II/base/logstream.h>

#include <vector>
#include <fstream>
#include <iomanip>

void f(const std::vector<int> &v)
{
  const VectorSlice<const std::vector<int> >
  s = make_slice(v,2,3);

  for (unsigned int i=0; i<s.size(); ++i)
    deallog << '\t' << s[i];
  deallog << std::endl;
}


int main()
{
  deal_II_exceptions::disable_abort_on_exception();

  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  std::vector<int> v(7);

  for (unsigned int i=0; i<v.size(); ++i)
    v[i] = i;

  VectorSlice<std::vector<int> > s(v, 3, 4);

  for (unsigned int i=0; i<s.size(); ++i)
    s[i] = i;

  for (unsigned int i=0; i<v.size(); ++i)
    deallog << '\t' << v[i];
  deallog << std::endl;

  f(v);

  try
    {
      make_slice(v, 3, 5);
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }
}
