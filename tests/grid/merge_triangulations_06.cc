// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 by the deal.II authors
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
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"


// Test the version of GridTools::merge_triangulations that takes a
// std::initializer_list.

void
test_2d()
{
  constexpr int      dim = 2;
  Triangulation<dim> tria_1;
  {
    GridGenerator::hyper_cube(tria_1);
    Tensor<1, dim> shift_vector;
    shift_vector[0] = -1.;
    shift_vector[1] = -1.;
    GridTools::shift(shift_vector, tria_1);
  }

  Triangulation<dim> tria_2;
  {
    GridGenerator::hyper_cube(tria_2);
    Tensor<1, dim> shift_vector;
    shift_vector[1] = -1.;
    GridTools::shift(shift_vector, tria_2);
  }

  Triangulation<dim> tria_3;
  {
    GridGenerator::hyper_cube(tria_3);
    Tensor<1, dim> shift_vector;
    shift_vector[0] = -1.;
    GridTools::shift(shift_vector, tria_3);
  }

  Triangulation<dim> tria_4;
  {
    GridGenerator::hyper_cube(tria_4);
  }

  // now merge triangulations
  Triangulation<dim> result;
  GridGenerator::merge_triangulations({&tria_1, &tria_2, &tria_3, &tria_4},
                                      result);

  GridOut().write_gnuplot(result, deallog.get_file_stream());

  deallog << "     Total number of cells        = " << result.n_cells()
          << std::endl
          << "     Total number of vertices = " << result.n_used_vertices()
          << std::endl;
}

void
test_3d()
{
  constexpr int      dim = 3;
  Triangulation<dim> tria_1;
  {
    GridGenerator::hyper_cube(tria_1);
    Tensor<1, dim> shift_vector;
    shift_vector[0] = -1.;
    shift_vector[1] = -1.;
    shift_vector[2] = -1.;
    GridTools::shift(shift_vector, tria_1);
  }

  Triangulation<dim> tria_2;
  {
    GridGenerator::hyper_cube(tria_2);
    Tensor<1, dim> shift_vector;
    shift_vector[1] = -1.;
    shift_vector[2] = -1.;
    GridTools::shift(shift_vector, tria_2);
  }

  Triangulation<dim> tria_3;
  {
    GridGenerator::hyper_cube(tria_3);
    Tensor<1, dim> shift_vector;
    shift_vector[0] = -1.;
    shift_vector[2] = -1.;
    GridTools::shift(shift_vector, tria_3);
  }

  Triangulation<dim> tria_4;
  {
    Tensor<1, dim> shift_vector;
    shift_vector[2] = -1.;
    GridGenerator::hyper_cube(tria_4);
    GridTools::shift(shift_vector, tria_4);
  }

  Triangulation<dim> tria_5;
  {
    GridGenerator::hyper_cube(tria_5);
    Tensor<1, dim> shift_vector;
    shift_vector[0] = -1.;
    shift_vector[1] = -1.;
    GridTools::shift(shift_vector, tria_5);
  }

  Triangulation<dim> tria_6;
  {
    GridGenerator::hyper_cube(tria_6);
    Tensor<1, dim> shift_vector;
    shift_vector[1] = -1.;
    GridTools::shift(shift_vector, tria_6);
  }

  Triangulation<dim> tria_7;
  {
    GridGenerator::hyper_cube(tria_7);
    Tensor<1, dim> shift_vector;
    shift_vector[0] = -1.;
    GridTools::shift(shift_vector, tria_7);
  }

  Triangulation<dim> tria_8;
  {
    GridGenerator::hyper_cube(tria_8);
  }

  // now merge triangulations
  Triangulation<dim> result;
  GridGenerator::merge_triangulations(
    {&tria_1, &tria_2, &tria_3, &tria_4, &tria_5, &tria_6, &tria_7, &tria_8},
    result);

  GridOut().write_gnuplot(result, deallog.get_file_stream());

  deallog << "     Total number of cells        = " << result.n_cells()
          << std::endl
          << "     Total number of vertices = " << result.n_used_vertices()
          << std::endl;
}

int
main()
{
  initlog();

  test_2d();
  test_3d();

  return 0;
}
