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

template <int dim, typename stream_type>
void
test(int n, const Point<dim> &split_point, stream_type &deallog)
{
  QSplit<dim> quad(QSimplex<dim>(QIterated<dim>(QTrapez<1>(), n)), split_point);

  deallog << std::endl
          << "# dim = " << dim << ", quad size = " << quad.size()
          << ", computed area = "
          << std::accumulate(quad.get_weights().begin(),
                             quad.get_weights().end(),
                             0.0)
          << std::endl
          << std::endl;

  for (auto p : quad.get_points())
    deallog << p << std::endl;
}


int
main()
{
  initlog();
  test<1>(4, Point<1>(.3), deallog);
  test<2>(4, Point<2>(.3, .4), deallog);
  test<3>(4, Point<3>(.3, .4, .5), deallog);
}
