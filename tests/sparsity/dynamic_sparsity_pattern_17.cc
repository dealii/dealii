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



// check DynamicSparsityPattern::column_index.
// this test is based on dynamic_sparsity_pattern_04.cc

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

  for (unsigned int i = 0; i < N; ++i)
    for (unsigned int ind = 0; ind < csp.row_length(i); ++ind)
      {
        const auto j   = csp.column_number(i, ind);
        const auto val = csp.column_index(i, j);
        AssertThrow(val == ind,
                    ExcMessage(std::to_string(val) +
                               "!=" + std::to_string(ind)));
      }

  deallog << "Ok" << std::endl;
}



int
main()
{
  initlog();

  test();
  return 0;
}
