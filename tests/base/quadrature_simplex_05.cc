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

template <int dim, typename stream_type>
void
test(int n, const Point<dim> &split_point, stream_type &deallog)
{
  QSplit<dim> quad(QSimplex<dim>(QIterated<dim>(QTrapezoid<1>(), n)),
                   split_point);

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
