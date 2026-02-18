// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2008 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// check DynamicSparsityPattern with an IndexSet that stores
// a non-contiguous range

#include "sparsity_pattern_common.h"

int
main()
{
  initlog();
  deallog << std::setprecision(3) << std::fixed;

  test_index_set<DynamicSparsityPattern>(false);
}
