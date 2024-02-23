// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <iostream>

#include "../tests.h"

// This test creates a parallelepiped (one element) and the same
// parallelepiped (one element); because the code that makes the two
// one-element grids is not the same, but should result in the same
// grid), this test is not entirely trivial.
//
// Here is the implementation in nd:
template <int dim>
void
check_nd_parallelepiped_by_comparison(bool log)
{
  // build corners for this particular dim that are known to give the
  // same output order as subdivided parallelepiped:
  Point<dim>(corners)[dim];

  switch (dim)
    {
      case 1:
        {
          corners[0] = Point<dim>(0.50);
          break;
        }

      case 2:
        {
          corners[0] = Point<dim>(0.50, 0.25);
          corners[1] = Point<dim>(0.25, 0.50);
          break;
        }

      case 3:
        {
          corners[0] = Point<dim>(0.25, 0.50, 0.50);
          corners[1] = Point<dim>(0.50, 0.25, 0.50);
          corners[2] = Point<dim>(0.50, 0.50, 0.25);
          break;
        }

      default:
        Assert(false, ExcInternalError());
    };

  Triangulation<dim> triangulation_parallelepiped;
  GridGenerator::parallelepiped(triangulation_parallelepiped, corners, false);

  Triangulation<dim> triangulation_subdivided_parallelepiped;
  GridGenerator::subdivided_parallelepiped(
    triangulation_subdivided_parallelepiped, 1, corners, false);

  if (log)
    {
      deallog << "check " << dim
              << "d parallelepiped (subdivided_parallelepiped):";
      if (GridTools::have_same_coarse_mesh(
            triangulation_parallelepiped,
            triangulation_subdivided_parallelepiped))
        deallog << "OK" << std::endl;

      else
        deallog
          << "not OK... coarse grids are different but they should be the same"
          << std::endl;
    }
}

int
main()
{
  initlog();
  // Check parallelepiped
  check_nd_parallelepiped_by_comparison<1>(true);
  check_nd_parallelepiped_by_comparison<2>(true);
  check_nd_parallelepiped_by_comparison<3>(true);
}
