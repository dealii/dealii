// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test IndexSet::get_index_vector()

#include <deal.II/base/index_set.h>

#include <stdlib.h>

#include "../tests.h"


void
test()
{
  IndexSet is1(100);

  // randomly add 90 elements to each
  // set, some of which may be
  // repetitions of previous ones
  for (unsigned int i = 0; i < 9 * is1.size() / 10; ++i)
    is1.add_index(Testing::rand() % is1.size());

  const std::vector<types::global_dof_index> indices = is1.get_index_vector();

  deallog << "Original index set: " << std::endl;
  is1.print(deallog);

  deallog << "List of indices: " << std::endl;
  for (unsigned int i = 0; i < indices.size(); ++i)
    deallog << indices[i] << ' ';
  deallog << std::endl;

  for (unsigned int i = 0; i < indices.size(); ++i)
    Assert(is1.index_within_set(indices[i]) == i, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
