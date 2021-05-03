// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2019 by the deal.II authors
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



// Test GridGenerator::eccentric_hyper_shell

#include <deal.II/base/tensor.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
test(std::ostream &out)
{
  Triangulation<dim> triangulation;

  const Point<dim> outer_center;
  Point<dim>       inner_center;
  inner_center[0]           = 0.2;
  const double inner_radius = 0.5;
  const double outer_radius = 1.0;

  GridGenerator::eccentric_hyper_shell(triangulation,
                                       inner_center,
                                       outer_center,
                                       inner_radius,
                                       outer_radius,
                                       dim == 2 ? 10 : 12);

  triangulation.refine_global(1);

  GridOut go;
  deallog.get_file_stream() << "Output mesh for dim = " << dim << std::endl
                            << std::endl;
  go.write_gnuplot(triangulation, out);
}


int
main()
{
  initlog();

  test<2>(deallog.get_file_stream());
  test<3>(deallog.get_file_stream());
}
