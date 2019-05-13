// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
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
