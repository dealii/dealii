// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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



// Test GridGenerator::hyper_ball_balanced

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
check_grid()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_ball_balanced(triangulation, Point<dim>(), 1.);
  deallog << "Number of cells: " << triangulation.n_cells() << std::endl;
  deallog << "Number of vertices: " << triangulation.n_vertices() << std::endl;
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


int
main()
{
  initlog();

  check_grid<2>();
  check_grid<3>();
}
