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

// construct a simplex quadrature, and check that we can get an affine
// transformation out of it.


#include <deal.II/base/quadrature_lib.h>

#include <numeric>

#include "../tests.h"

#include "simplex.h"

void
test(int n)
{
  const unsigned int dim = 2;

  QTrianglePolar quad(n);

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
  test(10);
}
