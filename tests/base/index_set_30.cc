// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// test IndexSet::index_within_set () for an empty index set

#include "../tests.h"

#include <deal.II/base/index_set.h>

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
