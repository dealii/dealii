// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test that iterating over a DynamicSparsityPattern with begin()/end() also
// works for indices exceeding 32 bit integers

#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include "../tests.h"


void
test()
{
  // take a bit more indices than we would have for 32 bit integers
  IndexSet set((1ULL << 32) + 4);

  // add range for the last six indices
  set.add_range(set.size() - 6, set.size());

  DynamicSparsityPattern sp;
  sp.reinit(set.size(), set.size(), set);

  for (const auto i : set)
    {
      sp.add(i, i);
      sp.add(i, i - 1);
      if (i + 1 < set.size())
        sp.add(i, i + 1);
    }

  deallog << "Sparsity pattern has " << sp.n_rows() << " rows" << std::endl;
  for (const auto row : set)
    {
      deallog << "Row " << row << " has the column indices: ";
      for (auto it = sp.begin(row); it != sp.end(row); ++it)
        {
          deallog << it->column() << ' ';
        }
      deallog << std::endl;
    }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
