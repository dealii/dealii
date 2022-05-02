// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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
