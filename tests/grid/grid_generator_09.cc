// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <iostream>

#include "../tests.h"

template <int dim>
void
check_rect1(unsigned int n, bool color, bool log)
{
  Point<dim>                left;
  Point<dim>                right;
  std::vector<unsigned int> subdivisions(dim);

  for (unsigned int d = 0; d < dim; ++d)
    {
      left[d]         = -1.;
      right[d]        = d + 2;
      subdivisions[d] = n * (d + 3);
    }
  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_rectangle(
    tria, subdivisions, left, right, color);

  GridOut grid_out;
  if (dim == 2)
    {
      if (log)
        grid_out.write_xfig(tria, deallog.get_file_stream());
      else
        grid_out.write_xfig(tria, std::cout);
    }
  else
    {
      if (log)
        grid_out.write_dx(tria, deallog.get_file_stream());
      else
        grid_out.write_dx(tria, std::cout);
    }
}


int
main()
{
  initlog();

  check_rect1<2>(1, true, true);
  check_rect1<2>(3, true, true);
  check_rect1<3>(1, true, true);
}
