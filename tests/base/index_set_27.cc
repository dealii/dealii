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

// Check that pop_back() for `IndexSet` works properly

#include <deal.II/base/index_set.h>

#include "../tests.h"

int
main()
{
  initlog();

  const unsigned int size = 100;

  IndexSet locally_owned(size);
  locally_owned.add_range(0, size);

  deallog << locally_owned.nth_index_in_set(5) << std::endl;
}
