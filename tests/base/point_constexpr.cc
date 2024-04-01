// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test for constexpr Point

#include <deal.II/base/point.h>

#include "../tests.h"


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  static_assert((Point<1>() / 1.).norm_square() == 0.);
  static_assert(Point<2>{0., 1.}.distance_square(Point<2>{}) == 1.);
  static_assert((-Point<3>{0., 0., 1.}).square() == 1.);

  deallog << "Ok" << std::endl;
}
