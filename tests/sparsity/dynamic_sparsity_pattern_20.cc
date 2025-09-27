// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check DynamicSparsityPattern::nonempty_columns()
// and nonempty_rows()

#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include "../tests.h"


void
test()
{
  const unsigned int     N = 100;
  DynamicSparsityPattern dsp(N, 2 * N);

  for (unsigned int i = 0; i < 5; ++i)
    for (unsigned int j = (i == 0 ? 0 : i - 1); j < i + 1; ++j)
      dsp.add(i, j);

  dsp.add(50, 50);
  dsp.add(51, 51);
  dsp.add(90, 52);

  dsp.add(N - 1, 2 * N - 1);

  const IndexSet cols = dsp.nonempty_cols();
  deallog << "Columns:" << std::endl;
  cols.print(deallog.get_file_stream());

  const IndexSet rows = dsp.nonempty_rows();
  deallog << "Rows:" << std::endl;
  rows.print(deallog.get_file_stream());
}



int
main()
{
  initlog();

  test();
  return 0;
}
