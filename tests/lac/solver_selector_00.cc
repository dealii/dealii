// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// I managed to break solver selector once by making some variables in the
// AdditionalData structures of the solvers const. This test simply
// instantiates that class, to make sure it still compiles

#include <deal.II/lac/solver_selector.h>

#include "../tests.h"

DEAL_II_NAMESPACE_OPEN
// instantiation here
template class SolverSelector<>;
DEAL_II_NAMESPACE_CLOSE

int
main()
{
  initlog();

  deallog << "OK" << std::endl;
}
