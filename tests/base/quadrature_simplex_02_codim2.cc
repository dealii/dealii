// ---------------------------------------------------------------------

// Copyright (C) 2022 by the deal.II authors

// This file is part of the deal.II library.

// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.

// ---------------------------------------------------------------------

// construct a simplex quadrature, and check that the affine transformation
// gives the correct results in the codimension one case

#include <deal.II/base/quadrature_lib.h>

#include <numeric>

#include "../tests.h"

// #include "simplex.h"

template <int dim, int spacedim>
void
test(const int degree, const std::array<Point<spacedim>, dim + 1> &vertices)
{
  QGaussSimplex<dim> quad(degree);

  deallog << "# dim = " << dim << std::endl;
  deallog << "# spacedim = " << spacedim << std::endl;

  auto quad2 = quad.compute_affine_transformation(vertices);

  for (auto p : quad2.get_points())
    deallog << p << std::endl;

  deallog << std::endl
          << "# Area: " << std::setprecision(15)
          << std::accumulate(quad2.get_weights().begin(),
                             quad2.get_weights().end(),
                             0.0)
          << std::endl
          << std::endl;
}


int
main()
{
  initlog();
  test<1, 3>(1, {{Point<3>(0., 1., 0.), Point<3>(0., 0., 1.)}});
}
