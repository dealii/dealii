// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check that the move constructor for `IndexSet` works properly

#include <deal.II/base/index_set.h>

#include "../tests.h"

int
main()
{
  initlog();

  IndexSet is1(100);
  is1.add_range(0, 10);
  is1.add_range(30, 40);

  deallog << is1.size() << ", " << is1.n_elements() << std::endl;

  // Test that move construction works correctly and that the moved object is
  // restored to the default state
  IndexSet is2 = std::move(is1);

  deallog << is2.size() << ", " << is2.n_elements() << std::endl;
  deallog << is1.size() << ", " << is1.n_elements() << std::endl;

  // Test that re-initializing the moved variable works
  is1.set_size(200);
  is1.add_range(90, 110);
  is1.add_range(130, 140);
  is1.add_range(145, 150);

  deallog << is1.size() << ", " << is1.n_elements() << std::endl;

  // Test that move assignment works correctly and that the moved object is
  // restored to the default state
  is2 = std::move(is1);
  deallog << is2.size() << ", " << is2.n_elements() << std::endl;
  deallog << is1.size() << ", " << is1.n_elements() << std::endl;
}
