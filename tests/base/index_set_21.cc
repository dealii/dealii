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


// test IndexSet::clear

#include <deal.II/base/index_set.h>

#include "../tests.h"


void
test()
{
  IndexSet index_set(20);
  index_set.add_range(2, 4);
  index_set.add_range(12, 18);
  index_set.add_index(6);
  index_set.add_index(8);
  index_set.add_index(14);
  index_set.add_index(16);

  // clear the IndexSet and then set elements
  // again
  index_set.clear();

  index_set.add_range(2, 4);
  index_set.add_range(12, 18);
  index_set.add_index(6);
  index_set.add_index(8);
  index_set.add_index(14);
  index_set.add_index(16);

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
