// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2018 by the deal.II authors
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


// test indexing in IndexSet variables for a non-contiguous range

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

  index_set.add_index(9);

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
