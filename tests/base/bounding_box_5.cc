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

// test BoundingBox::extend

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

  box.extend(.5);

  deallog << "New points: " << box.get_boundary_points().first << ", "
          << box.get_boundary_points().second << std::endl;


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
