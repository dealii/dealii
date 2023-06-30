// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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


// Test that IndexSet::add_index runs in linear complexity when adding the
// same indices many times

#include <deal.II/base/index_set.h>

#include "../tests.h"


void
test()
{
  IndexSet set(20);

  // Choose a large number of entries, which would time-out for quadratic
  // complexity
  for (unsigned int t = 0; t < 200000; ++t)
    for (unsigned int i = 4; i < 15; ++i)
      set.add_index(i);

  set.compress();
  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
