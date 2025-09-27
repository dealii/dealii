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

// test the implementation of the Utilities::mutable_bind function with
// std::bind objects


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
