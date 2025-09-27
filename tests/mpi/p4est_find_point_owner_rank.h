// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/mpi.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

// STL
#include <algorithm> // std::adjacent_find

/*
 * Test if the rank of point owners returned by all MPI processes is the same.
 */
template <int dim>
void
test(const std::vector<Point<dim>> &points)
{
  deallog << " Point search on subdivided hyper cube..." << std::endl;

  parallel::distributed::Triangulation<dim> triangulation(
    MPI_COMM_WORLD,
    typename Triangulation<dim>::MeshSmoothing(
      Triangulation<dim>::smoothing_on_refinement |
      Triangulation<dim>::smoothing_on_coarsening),
    typename parallel::distributed::Triangulation<dim>::Settings(
      parallel::distributed::Triangulation<
        dim>::communicate_vertices_to_p4est));
  GridGenerator::subdivided_hyper_cube(triangulation, 2);
  triangulation.refine_global(3);

  deallog << "   Number of MPI processes = "
          << Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) << std::endl;
  deallog << "   Number of this MPI processes = "
          << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) << std::endl;
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
  std::vector<types::subdomain_id> point_owner_ranks =
    triangulation.find_point_owner_rank(points);

  // Gather all point_owner_ranks found from all MPI ranks and compare the
  // results on root_process. They must all be equal.
  const unsigned int                            root_process{0};
  std::vector<std::vector<types::subdomain_id>> all_ranks_found_on_each_process{
    Utilities::MPI::gather(MPI_COMM_WORLD, point_owner_ranks, root_process)};

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == root_process)
    {
      if (std::adjacent_find(all_ranks_found_on_each_process.begin(),
                             all_ranks_found_on_each_process.end(),
                             std::not_equal_to<>()) ==
          all_ranks_found_on_each_process.end())
        {
          deallog << std::endl
                  << "   All point owner ranks found on all processes equal "
                     "each other. TEST PASSED."
                  << std::endl
                  << std::endl;
        }
      else
        {
          deallog
            << std::endl
            << "   Some point owner ranks found on all processes do not equal "
               "each other. Check output and/or the file "
               "tests/mpi/p4est_find_point_owner_rank.h. TEST FAILED."
            << std::endl
            << std::endl;
        }
    }
  ////////////////////////////////////////////////////////////

  deallog << "   >>> Reached end of test <<<" << std::endl;
}


/*
 * Vector input version
 */
template <int dim>
void
check_error_on_invalid_mesh(const Point<dim> &point)
{
  deallog << " Point search on hyper shell (manifold attached)..." << std::endl;

  parallel::distributed::Triangulation<dim> triangulation(
    MPI_COMM_WORLD,
    typename Triangulation<dim>::MeshSmoothing(
      Triangulation<dim>::smoothing_on_refinement |
      Triangulation<dim>::smoothing_on_coarsening),
    typename parallel::distributed::Triangulation<dim>::Settings(
      parallel::distributed::Triangulation<
        dim>::communicate_vertices_to_p4est));

  GridGenerator::hyper_shell(triangulation,
                             /* center is origin */ Point<dim>(),
                             /* 	inner_radius */ 1,
                             /* 	outer_radius */ 2);
  triangulation.refine_global(3);

  deallog << "   Number of MPI processes = "
          << Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) << std::endl;
  deallog << "   Number of this MPI processes = "
          << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) << std::endl;
  deallog << "   Number of cells = " << triangulation.n_global_active_cells()
          << std::endl;
  deallog << "   Number of levels = " << triangulation.n_levels() << std::endl;
  deallog << "   Number of global levels = " << triangulation.n_global_levels()
          << std::endl;
  const unsigned int checksum = triangulation.get_checksum();
  deallog << "   Triangulation checksum = " << checksum << std::endl;
  deallog << "   point = " << point << std::endl;

  deallog << std::endl;

  ////////////////////////////////////////////////////////////
  // This should result in an error.
  types::subdomain_id point_owner_rank =
    triangulation.find_point_owner_rank(point);
  ////////////////////////////////////////////////////////////
}
