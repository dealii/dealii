// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check DynamicSparsityPattern::copy operator with offdiagonals

#include "sparsity_pattern_common.h"

int
main()
{
  initlog();
  deallog << std::setprecision(3) << std::fixed;

  copy_with_offdiagonals_1<DynamicSparsityPattern>();
}
