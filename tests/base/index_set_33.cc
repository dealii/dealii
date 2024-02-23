// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
