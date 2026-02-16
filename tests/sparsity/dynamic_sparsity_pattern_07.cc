// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2004 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// check DynamicSparsityPattern::exists. since we create quite some
// output here, choose smaller number of rows and entries than in the other
// tests

#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include "../tests.h"


void
test()
{
  const unsigned int     N = 100;
  DynamicSparsityPattern csp(N, N);
  for (unsigned int i = 0; i < N; ++i)
    for (unsigned int j = 0; j < 10; ++j)
      csp.add(i, (i + (i + 1) * (j * j + i)) % N);
  csp.symmetrize();

  for (unsigned int i = 0; i < N; ++i)
    {
      deallog << i << ' ';
      for (unsigned int j = 0; j < N; ++j)
        deallog << (csp.exists(i, j) ? 'x' : ' ');
      deallog << std::endl;
    }
}



int
main()
{
  initlog();

  test();
  return 0;
}
