// ------------------------------------------------------------------------

// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors

// This file is part of the deal.II library.

// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.

// ------------------------------------------------------------------------

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
