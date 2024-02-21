// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test for Point::distance_square()

#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"


template <int dim>
void
check()
{
  Point<dim> p1, p2;
  for (unsigned int i = 0; i < dim; ++i)
    {
      p1[i] = 10.0 + 0.12345 * i;
      p1[i] = 0.5 + 0.6789 * i;
    }

  const double d  = p1.distance(p2);
  const double d2 = p1.distance_square(p2);

  AssertThrow(std::abs(d - std::sqrt(d2)) < 1e-10, ExcInternalError());

  deallog << "Ok" << std::endl;
}

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  check<1>();
  check<2>();
  check<3>();
}
