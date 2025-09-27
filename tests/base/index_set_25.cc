// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check that pop_front() for `IndexSet` works properly

#include <deal.II/base/index_set.h>

#include "../tests.h"

int
main()
{
  initlog();

  IndexSet is1(10);
  is1.add_range(0, 2);
  is1.add_range(5, 8);

  deallog << is1.n_elements() << ", ";
  deallog << is1.pop_front() << std::endl;

  deallog << is1.n_elements() << ", ";
  deallog << is1.pop_front() << std::endl;

  deallog << is1.n_elements() << ", ";
  deallog << is1.pop_front() << std::endl;


  is1.add_index(1);

  deallog << is1.n_elements() << ", ";
  deallog << is1.pop_front() << std::endl;

  deallog << is1.n_elements() << ", ";
  deallog << is1.pop_front() << std::endl;
}
