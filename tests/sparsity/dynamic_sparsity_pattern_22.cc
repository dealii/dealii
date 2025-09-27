// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check DynamicSparsityPattern::get_view()

#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include "../tests.h"


void
test()
{
  const unsigned int     N = 10;
  DynamicSparsityPattern dsp(N, N);

  for (unsigned int i = 0; i < N; ++i)
    if (i != 2 && i != 4)
      for (unsigned int j = (i == 0 ? 0 : i - 1); j < (i == N - 1 ? N : i + 1);
           ++j)
        dsp.add(i, j);

  IndexSet rows(N);
  rows.add_range(1, 5);

  DynamicSparsityPattern view = dsp.get_view(rows);

  deallog << "Initial sparsity:" << std::endl;
  dsp.print(deallog.get_file_stream());
  deallog << "Rows:" << std::endl;
  rows.print(deallog.get_file_stream());
  deallog << "View:" << std::endl;
  view.print(deallog.get_file_stream());
}



int
main()
{
  initlog();

  test();
  return 0;
}
