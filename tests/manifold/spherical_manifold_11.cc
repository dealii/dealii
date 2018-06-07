// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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

// test tangent vectors for nearly linear dependent positions

#include <deal.II/base/utilities.h>

#include <deal.II/grid/manifold_lib.h>

#include "../tests.h"


int
main()
{
  initlog();
  deallog << std::setprecision(10);

  const SphericalManifold<2> manifold_2;
  const SphericalManifold<3> manifold_3;

  // get tangent vectors for positions that are on a line and almost on a line
  // as seen from the center
  deallog
    << manifold_2.get_tangent_vector(Point<2>(0.5, 0), Point<2>(1.0, 0))
    << std::endl
    << manifold_2.get_tangent_vector(Point<2>(0.5, 0), Point<2>(1.0, 0.01))
    << std::endl
    << manifold_2.get_tangent_vector(Point<2>(0.5, 0), Point<2>(1.0, 1e-6))
    << std::endl
    << manifold_2.get_tangent_vector(Point<2>(0.5, 0), Point<2>(1.0, -1e-6))
    << std::endl
    << manifold_3.get_tangent_vector(Point<3>(0, 0, -0.5), Point<3>(0, 0, -1.0))
    << std::endl
    << manifold_3.get_tangent_vector(Point<3>(0, 0, -0.5),
                                     Point<3>(0, 1e-6, -1.0))
    << std::endl
    << manifold_3.get_tangent_vector(Point<3>(0, 0, 1.0),
                                     Point<3>(0, -1e-6, 0.5))
    << std::endl
    << manifold_3.get_tangent_vector(Point<3>(0.2, 0.3, 0.4),
                                     Point<3>(0.4, 0.55, 0.8))
    << std::endl
    << manifold_3.get_tangent_vector(
         Point<3>(std::sqrt(0.5), std::sqrt(0.5), 0.),
         Point<3>(-std::sqrt(0.5), std::sqrt(0.5), 0.))
    << std::endl
    << manifold_3.get_tangent_vector(
         Point<3>(std::sqrt(1. / 3.), std::sqrt(1. / 3.), std::sqrt(1. / 3.)),
         Point<3>(-std::sqrt(1. / 3.), std::sqrt(1. / 3.), std::sqrt(1. / 3.)))
    << std::endl
    << manifold_3.get_tangent_vector(
         Point<3>(std::sqrt(1. / 3.), std::sqrt(1. / 3.), std::sqrt(1. / 3.)),
         Point<3>(-std::sqrt(1. / 3.), -std::sqrt(1. / 3.), std::sqrt(1. / 3.)))
    << std::endl
    << manifold_3.get_tangent_vector(
         Point<3>(std::sqrt(1. / 3.), std::sqrt(1. / 3.), std::sqrt(1. / 3.)),
         Point<3>(std::sqrt(1. / 3.), -std::sqrt(1. / 3.), std::sqrt(1. / 3.)))
    << std::endl;
}
