// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// construct a simplex quadrature, and check that we can get an affine
// transformation out of it.

#include <deal.II/base/quadrature_lib.h>

#include <numeric>

#include "../tests.h"

#include "simplex.h"

template <int dim>
void
test(int n)
{
  QSimplex<dim> quad(QIterated<dim>(QTrapezoid<1>(), n));

  deallog << "# dim = " << dim << std::endl;

  for (auto p : quad.get_points())
    deallog << p << std::endl;

  deallog << std::endl
          << "# Area: "
          << std::accumulate(quad.get_weights().begin(),
                             quad.get_weights().end(),
                             0.0)
          << std::endl
          << std::endl;

  auto quad2 = quad.compute_affine_transformation(get_simplex<dim>());

  for (auto p : quad2.get_points())
    deallog << p << std::endl;


  deallog << std::endl
          << "# Area 2: "
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

  test<1>(10);
  test<2>(10);
  test<3>(10);
}
