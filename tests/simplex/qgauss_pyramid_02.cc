// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Tests QGaussPyramid by comparing it to the quadrautre points and weights
// of the old implementation.

#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include "../tests.h"

template <int dim>
void
test(const unsigned int n_points_1D)
{
  deallog << "N points 1D: " << n_points_1D << std::endl;

  QGaussPyramid<dim> quad(n_points_1D);
  for (unsigned int i = 0; i < quad.size(); ++i)
    deallog << quad.point(i) << " " << quad.weight(i) << std::endl;

  deallog << std::endl;
}


int
main()
{
  initlog();

  test<3>(1);
  test<3>(2);
}
