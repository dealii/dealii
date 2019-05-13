// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



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
