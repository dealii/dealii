// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// apply SparsityTools::reorder_hierarchical to a graph that consists
// of two non-connected parts.

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>

#include "../tests.h"


int
main()
{
  initlog();
  deallog << std::setprecision(3) << std::fixed;

  DynamicSparsityPattern dsp(8, 8);
  for (unsigned int i = 0; i < 8; ++i)
    dsp.add(i, i);

  // create a graph with components 0,2 and 1,3 that are disconnected
  dsp.add(0, 2);
  dsp.add(2, 0);

  dsp.add(1, 3);
  dsp.add(3, 1);

  // couple all indices between 3 and 7 except 4
  for (unsigned int i = 3; i < 7; ++i)
    for (unsigned int j = 3; j < 7; ++j)
      if (i != 4 && j != 4)
        dsp.add(i, j);

  // indices 4,7 couple
  dsp.add(4, 7);
  dsp.add(7, 4);

  // now find permutation
  std::vector<types::global_dof_index> permutation(8);
  SparsityTools::reorder_hierarchical(dsp, permutation);

  for (unsigned int i = 0; i < permutation.size(); ++i)
    deallog << permutation[i] << std::endl;
}
