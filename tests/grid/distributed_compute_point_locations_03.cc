// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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

// Test GridTools::compute_point_locations, in particular checks if the
// maps which are sent/received are correct

#include <deal.II/base/function_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



template <int dim>
void
test_distributed_cpt(unsigned int ref_cube)
{
  MPI_Comm     mpi_communicator = MPI_COMM_WORLD;
  unsigned int n_procs = Utilities::MPI::n_mpi_processes(mpi_communicator);
  unsigned int my_rank = Utilities::MPI::this_mpi_process(mpi_communicator);

  deallog << "Testing for dim = " << dim << " on " << n_procs << " processes"
          << std::endl;
  deallog << "Cube refinements: " << ref_cube << std::endl;

  // Creeating the cube on which to run distributed cpt loc
  parallel::distributed::Triangulation<dim> cube_d(mpi_communicator);
  GridGenerator::hyper_cube(cube_d);
  cube_d.refine_global(ref_cube);

  // We shall use the points from a shared grid so that each index is known
  std::vector<Point<dim>> test_points;
  Triangulation<dim>      cube;
  GridGenerator::hyper_cube(cube);
  cube.refine_global(ref_cube);
  for (auto &cell : cube.active_cell_iterators())
    test_points.emplace_back(cell->center());

  deallog << " Testing on " << test_points.size() << " points" << std::endl;

  // Computing bounding boxes describing the locally owned part of the mesh
  IteratorFilters::LocallyOwnedCell locally_owned_cell_predicate;
  std::vector<BoundingBox<dim>>     local_bbox =
    GridTools::compute_mesh_predicate_bounding_box(
      cube_d,
      std::function<bool(
        const typename Triangulation<dim>::active_cell_iterator &)>(
        locally_owned_cell_predicate),
      1,
      false,
      4);

  // Obtaining the global mesh description through an all to all communication
  std::vector<std::vector<BoundingBox<dim>>> global_bboxes;
  global_bboxes = Utilities::MPI::all_gather(mpi_communicator, local_bbox);

  // Initializing the cache
  GridTools::Cache<dim, dim> cache_d(cube_d);
  auto                       output_tuple =
    distributed_compute_point_locations(cache_d, test_points, global_bboxes);
  const auto &maps   = std::get<2>(output_tuple);
  const auto &points = std::get<3>(output_tuple);
  const auto &ranks  = std::get<4>(output_tuple);

  // Testing the results: if the map in maps is correct test_points[maps[i][j]]
  // == points[i][j]
  bool test_passed = true;

  for (unsigned int i = 0; i < points.size(); ++i)
    {
      for (unsigned int j = 0; j < maps[i].size(); ++j)
        if ((test_points[maps[i][j]] - points[i][j]).norm() > 1e-10)
          {
            deallog << " Error in cell " << i << " with position " << j
                    << std::endl;
            deallog << " Received map was: " << maps[i][j] << std::endl;
            deallog << " From rank: " << ranks[i][j] << std::endl;
            test_passed = false;
          }
    }

  if (test_passed)
    deallog << "Test passed" << std::endl;
  else
    deallog << "Test FAILED" << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog << "Deal.II GridTools::distributed_compute_point_locations"
          << std::endl;
  deallog << "Test 3: maps values" << std::endl;
  deallog << "2D tests:" << std::endl;
  test_distributed_cpt<2>(3);
  deallog << "3D tests" << std::endl;
  test_distributed_cpt<3>(3);
}
