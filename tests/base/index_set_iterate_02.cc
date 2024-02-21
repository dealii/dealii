// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test IndexSet::at()

#include <deal.II/base/index_set.h>

#include "../tests.h"

void
test(IndexSet &index_set, unsigned int n)
{
  deallog << "n=" << n;

  IndexSet::ElementIterator it = index_set.at(n);

  deallog << " end?" << (it == index_set.end());
  if (it != index_set.end())
    {
      deallog << " value=" << *it;
    }

  deallog << std::endl;
}

void
test()
{
  IndexSet index_set(20);
  index_set.add_range(2, 5);
  index_set.add_index(9);

  index_set.print(deallog);

  test(index_set, 0);
  test(index_set, 2);
  test(index_set, 3);
  test(index_set, 4);
  test(index_set, 7);
  test(index_set, 9);
  test(index_set, 15);

  IndexSet empty(42);
  empty.print(deallog);
  test(empty, 6);
}



int
main()
{
  initlog();

  test();
}
