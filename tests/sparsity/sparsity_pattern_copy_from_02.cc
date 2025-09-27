// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// SparsityPattern::copy_from used to be of quadratic complexity for
// DynamicSparsityPattern arguments with many empty rows. Create a big
// sparsity pattern of 10m elements and only put an entry into the last
// row. The old code would time out.

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <list>
#include <set>

#include "../tests.h"


int
main()
{
  initlog();
  deallog << std::setprecision(3) << std::fixed;

  DynamicSparsityPattern csp(10000000, 100);
  csp.add(9999999, 1);

  SparsityPattern sp;
  sp.copy_from(csp);

  deallog << "OK" << std::endl;
}
