// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2017 by the deal.II authors
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



// check DynamicSparsityPattern::bandwidth

#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include "../tests.h"


void
test()
{
  // set up a sparsity pattern. since
  // DynamicSparsityPatterns are most
  // often used for 3d, use a rather large
  // number of entries per row
  const unsigned int     N = 1000;
  DynamicSparsityPattern csp(N, N);
  for (unsigned int i = 0; i < N; ++i)
    for (unsigned int j = 0; j < 40; ++j)
      csp.add(i, (i + (i + 1) * (j * j + i)) % N);

  deallog << csp.bandwidth() << std::endl;
}



int
main()
{
  initlog();

  test();
  return 0;
}
