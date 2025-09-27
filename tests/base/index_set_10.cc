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


// test indexing in IndexSet variables for a contiguous range

#include <deal.II/base/index_set.h>

#include "../tests.h"


void
test()
{
  IndexSet index_set(20);

  index_set.add_index(2);
  index_set.add_index(3);
  index_set.add_index(4);

  index_set.add_index(6);
  index_set.add_index(7);

  index_set.compress();

  index_set.add_index(5);

  for (unsigned int i = 0; i < index_set.n_elements(); ++i)
    {
      deallog << index_set.nth_index_in_set(i) << std::endl;
      AssertThrow(index_set.index_within_set(index_set.nth_index_in_set(i)) ==
                    i,
                  ExcInternalError());
    }
  deallog << "OK" << std::endl;

  for (unsigned int i = 0; i < index_set.size(); ++i)
    if (index_set.is_element(i))
      deallog << i << ' ' << index_set.index_within_set(i) << std::endl;
}



int
main()
{
  initlog();

  test();
}
