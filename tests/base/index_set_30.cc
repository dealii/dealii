// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// test IndexSet::index_within_set () for an empty index set

#include <deal.II/base/index_set.h>

#include "../tests.h"


void
test()
{
  IndexSet index_set(20);

  index_set.compress();

  AssertThrow(index_set.index_within_set(2) == numbers::invalid_dof_index,
              ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
