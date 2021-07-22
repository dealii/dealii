// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2021 by the deal.II authors
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

#pragma once

#include <deal.II/base/mpi.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

/*
 * Single point version
 */
template <int dim>
void
test(const Point<dim> &point)
{
  deallog << " Point search on subdivided hyper cube..." << std::endl;

  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  GridGenerator::subdivided_hyper_cube(triangulation, 2);
  triangulation.refine_global(3);

  deallog << "   Number of cells = " << triangulation.n_global_active_cells()
          << std::endl;
  deallog << "   Number of levels = " << triangulation.n_levels() << std::endl;
  deallog << "   Number of global levels = " << triangulation.n_global_levels()
          << std::endl;
  const unsigned int checksum = triangulation.get_checksum();
  deallog << "   Triangulation checksum = " << checksum << std::endl;
  deallog << "   point = " << point << std::endl;

  ////////////////////////////////////////////////////////////
  // test stuff
  unsigned int point_owner_rank = triangulation.find_point_owner_rank(point);

  deallog << "   rank of point owner = " << point_owner_rank << std::endl;
  ////////////////////////////////////////////////////////////

  deallog << "   >>> Reached end of test <<<" << std::endl;
}



/*
 * Vector input version
 */
template <int dim>
void
test(const std::vector<Point<dim>> &points)
{
  deallog << " Point search on subdivided hyper cube..." << std::endl;

  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  GridGenerator::subdivided_hyper_cube(triangulation, 2);
  triangulation.refine_global(3);

  deallog << "   Number of cells = " << triangulation.n_global_active_cells()
          << std::endl;
  deallog << "   Number of levels = " << triangulation.n_levels() << std::endl;
  deallog << "   Number of global levels = " << triangulation.n_global_levels()
          << std::endl;
  const unsigned int checksum = triangulation.get_checksum();
  deallog << "   Triangulation checksum = " << checksum << std::endl;

  deallog << "   points = ";
  for (const auto &point : points)
    deallog << point << "    ";

  deallog << std::endl;

  ////////////////////////////////////////////////////////////
  // test stuff
  std::vector<unsigned int> point_owner_ranks =
    triangulation.find_point_owner_rank(points);

  deallog << "   rank of point owners = ";
  for (const auto &point_owner_rank : point_owner_ranks)
    deallog << point_owner_rank << "    ";

  deallog << std::endl;
  ////////////////////////////////////////////////////////////

  deallog << "   >>> Reached end of test <<<" << std::endl;
}
