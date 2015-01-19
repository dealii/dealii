// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2014 by the deal.II authors
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



// check SparsityPattern::copy_from(CompressedSetSparsityPattern)

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/compressed_set_sparsity_pattern.h>
#include <iomanip>
#include <fstream>


void test ()
{
  const unsigned int N = 100;
  CompressedSetSparsityPattern csp (N,N);
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<10; ++j)
      csp.add (i, (i+(i+1)*(j*j+i))%N);

  SparsityPattern sp;
  sp.copy_from (csp);
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<sp.row_length (i); ++j)
      deallog << i << ' ' << j << ' ' << sp.column_number (i,j)
              << std::endl;
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
