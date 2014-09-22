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



// check CompressedSparsityPattern::row_length

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <iomanip>
#include <fstream>


void test ()
{
  // set up a sparsity pattern. since
  // CompressedSparsityPatterns are most
  // often used for 3d, use a rather large
  // number of entries per row
  const unsigned int N = 1000;
  CompressedSparsityPattern csp (N,N);
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<40; ++j)
      csp.add (i, (i+(i+1)*(j*j+i))%N);

  for (unsigned int i=0; i<N; ++i)
    {
      std::vector<bool> xx(1000, false);
      for (unsigned int j=0; j<40; ++j)
        xx[(i+(i+1)*(j*j+i))%N] = true;

      Assert (static_cast<unsigned int>(std::count (xx.begin(), xx.end(),
                                                    true))
              == csp.row_length(i),
              ExcInternalError());

      deallog << i << ' ' << csp.row_length (i) << std::endl;
    }
}



int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
  return 0;
}
