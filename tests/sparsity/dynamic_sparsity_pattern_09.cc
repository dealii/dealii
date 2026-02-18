// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2004 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// check DynamicSparsityPattern::empty

#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include "../tests.h"


void
test()
{
  const unsigned int     N = 1000;
  DynamicSparsityPattern csp;
  AssertThrow(csp.empty() == true, ExcInternalError());

  csp.reinit(N, N);
  AssertThrow(csp.empty() == false, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
  return 0;
}
