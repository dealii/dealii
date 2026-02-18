// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2010 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// test for Point::operator()

#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"


template <int dim>
void
check()
{
  Point<dim> p;
  for (unsigned int i = 0; i < dim; ++i)
    p[i] = i;

  for (unsigned int i = 0; i < dim; ++i)
    deallog << p[i] << ' ';
  deallog << std::endl;
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
