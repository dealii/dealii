// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2017 by the deal.II authors
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


// test parallel::accumulate_from_subranges

#include <deal.II/base/parallel.h>

#include "../tests.h"


int
sum(const int begin, const int end)
{
  int s = 0;
  for (int i = begin; i < end; ++i)
    s += i;
  return s;
}


int
main()
{
  initlog();

  const int N = 10000;
  const int s = parallel::accumulate_from_subranges<int>(&sum, 0, N, 10);

  DEAL_II_Assert(s == N * (N - 1) / 2, ExcInternalError());

  deallog << s << std::endl;
}
