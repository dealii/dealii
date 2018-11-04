// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


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
