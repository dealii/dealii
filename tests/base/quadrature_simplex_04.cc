// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// construct an anisotropic simplex quadrature, and check that we can
// get an affine transformation out of it.


#include <deal.II/base/quadrature_lib.h>

#include <numeric>

#include "../tests.h"

#include "simplex.h"

void
test(unsigned int n)
{
  const unsigned int dim = 2;

  {
    QSplit<2> quad(QTrianglePolar(n), Point<2>(.3, .4));

    for (auto p : quad.get_points())
      deallog << p << std::endl;

    deallog << std::endl
            << "# Area: "
            << std::accumulate(quad.get_weights().begin(),
                               quad.get_weights().end(),
                               0.0)
            << std::endl
            << std::endl;

    if (quad.size() == n * n * 4)
      deallog << "# Size OK" << std::endl;
    else
      deallog << "# Size NOT OK" << std::endl;
  }
  {
    QSplit<2> quad(QTrianglePolar(n), Point<2>(0, .2));

    for (auto p : quad.get_points())
      deallog << p << std::endl;

    deallog << std::endl
            << "# Area: "
            << std::accumulate(quad.get_weights().begin(),
                               quad.get_weights().end(),
                               0.0)
            << std::endl
            << std::endl;
    if (quad.size() == n * n * 3)
      deallog << "# Size OK" << std::endl;
    else
      deallog << "# Size NOT OK" << std::endl;
  }
  {
    QSplit<2> quad(QTrianglePolar(n), Point<2>(1, 0));

    for (auto p : quad.get_points())
      deallog << p << std::endl;

    deallog << std::endl
            << "# Area: "
            << std::accumulate(quad.get_weights().begin(),
                               quad.get_weights().end(),
                               0.0)
            << std::endl
            << std::endl;
    if (quad.size() == n * n * 2)
      deallog << "# Size OK" << std::endl;
    else
      deallog << "# Size NOT OK" << std::endl;
  }
}


int
main()
{
  initlog();
  test(5);
}
