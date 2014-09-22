// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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


// test parallel::accumulate_from_subranges

#include "../tests.h"
#include <iomanip>
#include <fstream>

#include <deal.II/base/parallel.h>


int sum (const int begin,
         const int end)
{
  int s=0;
  for (int i=begin; i<end; ++i)
    s += i;
  return s;
}


int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  const int N = 10000;
  const int s = parallel::accumulate_from_subranges<int> (&sum, 0, N, 10);

  Assert (s == N*(N-1)/2, ExcInternalError());

  deallog << s << std::endl;
}
