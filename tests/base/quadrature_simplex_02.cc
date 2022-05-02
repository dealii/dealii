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
