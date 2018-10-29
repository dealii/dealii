// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2017 by the deal.II authors
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


// test IndexSet::n_elements()

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

  index_set.add_index(9);

  deallog << index_set.n_elements() << std::endl;
}



int
main()
{
  initlog();

  test();
}
