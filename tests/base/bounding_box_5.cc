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

// test BoundingBox::extend() and BoundingBox::create_extend()

#include <deal.II/base/bounding_box.h>

#include <iostream>
#include <vector>

#include "../tests.h"



template <int dim>
void
test()
{
  Point<dim> p0, p1;
  for (unsigned int d = 0; d < dim; ++d)
    p1[d] = 1;

  BoundingBox<dim> box({p0, p1});

  const auto box1 = box.create_extended(.5);

  deallog << "New points: " << box1.get_boundary_points().first << ", "
          << box1.get_boundary_points().second << std::endl;

  box.extend(.5);

  deallog << "New points: " << box.get_boundary_points().first << ", "
          << box.get_boundary_points().second << std::endl;

  const auto box2 = box.create_extended(-.8);

  deallog << "New points: " << box2.get_boundary_points().first << ", "
          << box2.get_boundary_points().second << std::endl;

  box.extend(-.8);

  deallog << "New points: " << box.get_boundary_points().first << ", "
          << box.get_boundary_points().second << std::endl;
}



int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();
}
