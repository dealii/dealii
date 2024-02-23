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



// investigate performance issues in DynamicSparsityPattern::begin(r) and
// end(r) for large sets with many empty rows.

#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include "../tests.h"


void
test(bool empty, bool large_gap)
{
  const int size     = 100000000;
  const int my_start = size / 3;
  IndexSet  owned(size);
  owned.add_range(my_start, my_start + 5);
  if (large_gap)
    owned.add_range(size - 1, size);
  DynamicSparsityPattern sp(size, 5, owned);
  if (!empty)
    sp.add(my_start + 1, 1);

  for (unsigned int i = my_start - 10; i < my_start + 10; ++i)
    for (DynamicSparsityPattern::iterator p = sp.begin(i); p != sp.end(i); ++p)
      deallog << p->row() << ' ' << p->column() << std::endl;

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test(false, false);
  test(true, false);
  test(false, true);
  test(true, true);
}
