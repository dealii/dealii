// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// test the implementation of the Utilities::mutable_bind function


#include <deal.II/base/mutable_bind.h>
#include <deal.II/base/point.h>

#include <tuple>

#include "../tests.h"

void
example_function(const Point<2> &p, const double &d, const unsigned int i = 3)
{
  deallog << "P: " << p << ", d: " << d << ", i: " << i << std::endl;
}

int
main()
{
  initlog();

  auto exp = Utilities::mutable_bind(example_function);
  exp();

  exp.set_arguments(std::make_tuple(Point<2>(.5, .5), 2.0, 3));
  exp();

  exp.parse_arguments("2.0 , 2.0 : 1.0 : 4");
  exp();

  exp.set_arguments(Point<2>(.3, .3), 2.0, 5);
  exp();

  auto exp1 =
    Utilities::mutable_bind(example_function, Point<2>(.8, .8), 4.0, 10);
  exp1();
}
