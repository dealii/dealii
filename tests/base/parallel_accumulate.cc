// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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

  Assert(s == N * (N - 1) / 2, ExcInternalError());

  deallog << s << std::endl;
}
