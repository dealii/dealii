// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2021 by the deal.II authors
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

// test the implementation of the Utilities::mutable_bind function with
// std::bind objects


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
  initlog(true);

  const Point<2> p(1, 2);

  Utilities::MutableBind<void, double, unsigned int> exp = {
    std::bind(example_function,
              std::cref(p),
              std::placeholders::_1,
              std::placeholders::_2),
    {}};

  exp();
  exp.set_arguments(3.0, 4);
  exp();
}
