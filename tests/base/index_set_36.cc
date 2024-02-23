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


// test IndexSet::get_view(IndexSet) with a mask that is identical to
// the one used in the _18 test just that it is given as an IndexSet
// rather than a single range.

#include <deal.II/base/index_set.h>

#include <stdlib.h>

#include "../tests.h"


void
test()
{
  IndexSet is1(100);

  // randomly add 90 elements to each
  // set, some of which may be
  // repetitions of previous ones
  for (unsigned int i = 0; i < 9 * is1.size() / 10; ++i)
    is1.add_index(Testing::rand() % is1.size());

  is1.compress();

  IndexSet mask(100);
  mask.add_range(20, 50);

  const IndexSet is2 = is1.get_view(mask);
  Assert(is2.size() == mask.n_elements(), ExcInternalError());

  deallog << "Original index set: " << std::endl;
  is1.print(deallog);
  deallog << "View of index set between 20 and 50: " << std::endl;
  is2.print(deallog);

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
