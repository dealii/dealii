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



// apply SparsityTools::reorder_Cuthill_McKee to a graph that consists
// of two or more non-connected parts. the reordering algorithm used
// to trip over that

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>

#include "../tests.h"


int
main()
{
  initlog();
  deallog << std::setprecision(3) << std::fixed;

  DynamicSparsityPattern dsp(4, 4);
  for (unsigned int i = 0; i < 4; ++i)
    dsp.add(i, i);

  // create a graph with components 0,2 and 1,3 that are disconnected
  dsp.add(0, 2);
  dsp.add(2, 0);

  dsp.add(1, 3);
  dsp.add(3, 1);

  // now find permutation
  std::vector<types::global_dof_index> permutation(4);
  SparsityTools::reorder_Cuthill_McKee(dsp, permutation);

  for (unsigned int i = 0; i < permutation.size(); ++i)
    deallog << permutation[i] << std::endl;
}
