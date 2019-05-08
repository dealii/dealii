// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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
