// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check that IncrementalFunction use correct number of components

#include <deal.II/base/incremental_function.h>

#include "../tests.h"

template <int dim>
void
check(const unsigned int n_comp = 3)
{
  Functions::ZeroFunction<dim>        func(n_comp);
  Functions::IncrementalFunction<dim> inc(func);

  AssertThrow(inc.n_components == n_comp, ExcInternalError());

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();

  check<2>();
}
