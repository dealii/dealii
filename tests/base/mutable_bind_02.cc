// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

// test the implementation of the Utilities::mutable_bind function. Make
// sure it works for std::function objects.

#include <deal.II/base/mutable_bind.h>
#include <deal.II/base/point.h>
#include <deal.II/base/std_cxx17/tuple.h>

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

  Point<2> p(.1, .1);

  std::function<void(const double, const unsigned int)> reduced_function =
    [&p](const double d, const unsigned int i) { example_function(p, d, i); };

  auto exp = Utilities::mutable_bind(reduced_function);
  exp();

  exp.set_arguments(std::make_tuple(2.0, 3));
  exp();

  exp.parse_arguments("1.0 : 4");
  exp();

  exp.set_arguments(2.0, 5);
  exp();

  auto exp1 = Utilities::mutable_bind(reduced_function, 4.0, 10);
  exp1();
}
