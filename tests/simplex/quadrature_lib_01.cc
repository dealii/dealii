// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test QGaussSimplex: output its quadrature points and weights.


#include <deal.II/base/quadrature_lib.h>

#include "../tests.h"


template <int dim>
void
test(const unsigned int n_points)
{
  QGaussSimplex<dim> quad(n_points);

  for (unsigned int q = 0; q < quad.size(); ++q)
    {
      deallog << quad.point(q) << ' ';
      deallog << quad.weight(q) << ' ';
      deallog << std::endl;
    }
}

int
main()
{
  initlog();

  {
    deallog.push("1d-1");
    test<1>(1 /*n_points*/);
    deallog.pop();
  }
  {
    deallog.push("1d-2");
    test<1>(2);
    deallog.pop();
  }
  {
    deallog.push("1d-3");
    test<1>(3);
    deallog.pop();
  }
  {
    deallog.push("2d-1");
    test<2>(1);
    deallog.pop();
  }
  {
    deallog.push("2d-3");
    test<2>(2);
    deallog.pop();
  }
  {
    deallog.push("2d-7");
    test<2>(3);
    deallog.pop();
  }

  {
    deallog.push("3d-1");
    test<3>(1);
    deallog.pop();
  }
  {
    deallog.push("3d-4");
    test<3>(2);
    deallog.pop();
  }
  {
    deallog.push("3d-10");
    test<3>(3);
    deallog.pop();
  }
}
