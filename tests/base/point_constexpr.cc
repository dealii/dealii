// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2010 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// test for constexpr Point

#include <deal.II/base/point.h>

#include "../tests.h"


template <int dim>
constexpr bool val =
  Point<dim>((Point<dim>::unit_vector(0) + Point<dim>::unit_vector(0)) -
             2. * Point<dim>::unit_vector(0))[0] == 0.;


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  static_assert(val<1>);
  static_assert(val<2>);
  static_assert(val<3>);

  static_assert((Point<1>() / 1.).norm_square() == 0.);
  static_assert(Point<2>{0., 1.}.distance_square(Point<2>{}) == 1.);
  static_assert((-Point<3>{0., 0., 1.}).square() == 1.);

  deallog << "Ok" << std::endl;
}
