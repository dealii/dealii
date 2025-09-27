// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



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
