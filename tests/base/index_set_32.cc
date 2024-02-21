// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test how IndexSet::add_indices deals with duplicate indices to be added

#include <deal.II/base/index_set.h>

#include "../tests.h"


void
test()
{
  IndexSet index_set(10);

  {
    const unsigned int array[] = {2, 2, 2, 3, 4, 6, 8, 8};
    index_set.add_indices((const unsigned int *)array,
                          array + sizeof(array) / sizeof(array[0]));
  }

  Assert(index_set.is_contiguous() == false, ExcInternalError());

  for (unsigned int i = 0; i < index_set.size(); ++i)
    deallog << i << ' ' << (index_set.is_element(i) ? "true" : "false")
            << std::endl;
}



int
main()
{
  initlog();

  test();
}
