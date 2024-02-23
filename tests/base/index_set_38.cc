// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test IndexSet::get_view(IndexSet) with a mask that is not just a
// single range.
//
// This test is for the specific example given in the documentation of
// the function we test.

#include <deal.II/base/index_set.h>

#include <stdlib.h>

#include "../tests.h"


void
test()
{
  IndexSet is1(100);
  is1.add_range(20, 40);
  is1.compress();

  IndexSet mask(100);
  for (unsigned int i = 0; i < is1.size(); ++i)
    if (i % 2 == 1)
      mask.add_index(i);
  mask.compress();

  const IndexSet is2 = is1.get_view(mask);
  Assert(is2.size() == mask.n_elements(), ExcInternalError());

  deallog << "Original index set: " << std::endl;
  is1.print(deallog);
  deallog << "Mask: " << std::endl;
  mask.print(deallog);
  deallog << "View of index set among the odd numbers: " << std::endl;
  is2.print(deallog);

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
