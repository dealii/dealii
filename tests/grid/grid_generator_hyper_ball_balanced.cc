// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test GridGenerator::hyper_ball_balanced

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
check_grid()
{
  // Check centered at the origin and at another point
  std::array<Point<dim>, 2> centers;
  centers[1][0]       = 0.6;
  centers[1][1]       = 0.5;
  centers[1][dim - 1] = 0.5;
  for (const Point<dim> &center : centers)
    {
      Triangulation<dim> triangulation;
      GridGenerator::hyper_ball_balanced(triangulation, center, 1.);
      deallog << "Number of cells: " << triangulation.n_cells() << std::endl;
      deallog << "Number of vertices: " << triangulation.n_vertices()
              << std::endl;
      triangulation.refine_global();

      GridOut            go;
      GridOutFlags::XFig xfig_flags;
      xfig_flags.fill_style = 25;

      go.set_flags(xfig_flags);
      GridOut::OutputFormat format = GridOut::xfig;
      if (dim == 3)
        format = GridOut::dx;

      go.write(triangulation, deallog.get_file_stream(), format);
    }
}


int
main()
{
  initlog();

  check_grid<2>();
  check_grid<3>();
}
