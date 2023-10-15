// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2018 by the deal.II authors
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


// test IndexSet::get_view(IndexSet) with a mask that is not just a
// single range.

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

  IndexSet mask(100);
  mask.add_range(20, 40);
  mask.add_range(60, 80);
  mask.compress();

  const IndexSet is2 = is1.get_view(mask);
  Assert(is2.size() == mask.n_elements(), ExcInternalError());

  deallog << "Original index set: " << std::endl;
  is1.print(deallog);
  deallog << "Mask: " << std::endl;
  mask.print(deallog);
  deallog << "View of index set between 20-40 and 60-80: " << std::endl;
  is2.print(deallog);

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
