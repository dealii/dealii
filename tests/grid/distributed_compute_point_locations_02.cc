// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test GridTools::distributed_compute_point_locations for the parallel case:
// Inside a distributed hypercube there's a shared sphere:
// call distributed point locations on the sphere's cells centers and check
// the result.

#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>

#include <deal.II/distributed/shared_tria.h>
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
test_compute_pt_loc(unsigned int ref_cube, unsigned int ref_sphere)
{
  MPI_Comm     mpi_communicator = MPI_COMM_WORLD;
  unsigned int n_procs = Utilities::MPI::n_mpi_processes(mpi_communicator);
  unsigned int my_rank = Utilities::MPI::this_mpi_process(mpi_communicator);

  deallog << "Testing for dim = " << dim << " on " << n_procs << " processes"
          << std::endl;
  deallog << "Cube refinements: " << ref_cube << std::endl;
  deallog << "Sphere refinements:" << ref_sphere << std::endl;

  // Initializing and refining meshes
  parallel::distributed::Triangulation<dim> cube(mpi_communicator);
  GridGenerator::hyper_cube(cube);
  cube.refine_global(ref_cube);

  parallel::shared::Triangulation<dim> sphere(mpi_communicator);
  Point<dim>                           sphere_center;
  // Defining center and radius
  for (unsigned int i = 0; i < dim; ++i)
    sphere_center[i] = 0.47 - i * 0.05;
  double radius = 0.4 - dim * 0.05;
  GridGenerator::hyper_ball(sphere, sphere_center, radius);
  static SphericalManifold<dim, dim> surface_description(sphere_center);
  sphere.set_manifold(0, surface_description);
  sphere.refine_global(ref_sphere);

  deallog << "Sphere center:" << sphere_center << std::endl;
  deallog << "Sphere radius:" << radius << std::endl;

  // Initializing the cache
  GridTools::Cache<dim, dim> cache(cube);

  // Centers of locally owned cells
  std::vector<Point<dim>> loc_owned_points;
  // Building by hand the output of distributed points (see the function's
  // description for more details, this is in fact compute point location
  // code with the addition of rank storing)
  std::vector<typename Triangulation<dim, dim>::active_cell_iterator>
                                         computed_cells;
  std::vector<std::vector<Point<dim>>>   computed_qpoints;
  std::vector<std::vector<Point<dim>>>   computed_points;
  std::vector<std::vector<unsigned int>> computed_ranks;

  unsigned int computed_pts = 0;
  for (auto &cell : sphere.active_cell_iterators())
    {
      // The points we consider are the cell centers
      auto center_pt = cell->center();
      // Store the point only if it is inside a locally owned sphere cell
      if (cell->subdomain_id() == my_rank)
        loc_owned_points.emplace_back(center_pt);
      // Find the cube cell where center pt lies
      auto my_pair = GridTools::find_active_cell_around_point(cache, center_pt);
      // If it is inside a locally owned cell it shall be returned
      // from distributed compute point locations
      if (my_pair.first.state() == IteratorState::valid &&
          my_pair.first->is_locally_owned())
        {
          computed_pts++;
          auto cells_it = std::find(computed_cells.begin(),
                                    computed_cells.end(),
                                    my_pair.first);

          if (cells_it == computed_cells.end())
            {
              // Cell not found: adding a new cell
              computed_cells.emplace_back(my_pair.first);
              computed_qpoints.emplace_back(1, my_pair.second);
              computed_points.emplace_back(1, center_pt);
              computed_ranks.emplace_back(1, cell->subdomain_id());
            }
          else
            {
              // Cell found: just adding the point index and qpoint to the
              // list
              unsigned int current_cell = cells_it - computed_cells.begin();
              computed_qpoints[current_cell].emplace_back(my_pair.second);
              computed_points[current_cell].emplace_back(center_pt);
              computed_ranks[current_cell].emplace_back(cell->subdomain_id());
            }
        }
    }

  // Computing bounding boxes describing the locally owned part of the mesh
  IteratorFilters::LocallyOwnedCell locally_owned_cell_predicate;
  std::vector<BoundingBox<dim>>     local_bbox =
    GridTools::compute_mesh_predicate_bounding_box(
      cache.get_triangulation(),
      std::function<bool(
        const typename Triangulation<dim>::active_cell_iterator &)>(
        locally_owned_cell_predicate),
      1,
      true,
      4);

  // Obtaining the global mesh description through an all to all communication
  std::vector<std::vector<BoundingBox<dim>>> global_bboxes;
  global_bboxes = Utilities::MPI::all_gather(mpi_communicator, local_bbox);

  // Using the distributed version of compute point location
  auto output_tuple =
    distributed_compute_point_locations(cache, loc_owned_points, global_bboxes);
  deallog << "Comparing results" << std::endl;
  const auto &output_cells   = std::get<0>(output_tuple);
  const auto &output_qpoints = std::get<1>(output_tuple);
  const auto &output_points  = std::get<3>(output_tuple);
  const auto &output_ranks   = std::get<4>(output_tuple);

  // Comparing the output with the previously computed computed result
  bool test_passed = true;
  if (output_cells.size() != computed_cells.size())
    {
      test_passed = false;
      deallog << "ERROR: non-matching number of cell found" << std::endl;
    }

  unsigned int output_computed_pts = 0;
  for (unsigned int c = 0; c < output_cells.size(); ++c)
    {
      output_computed_pts += output_points[c].size();
      const auto &cell = output_cells[c];
      auto        cell_it =
        std::find(computed_cells.begin(), computed_cells.end(), cell);
      if (cell_it == computed_cells.end())
        {
          deallog << "ERROR: active cell " << cell->active_cell_index()
                  << " not found" << std::endl;
          test_passed = false;
        }
      else
        {
          unsigned int c_cell = cell_it - computed_cells.begin();
          if (output_points[c].size() != computed_points[c_cell].size())
            {
              test_passed = false;
              deallog << "ERROR: non-matching number of points for cell "
                      << cell->active_cell_index() << std::endl;
              deallog << "Distributed compute point location output:"
                      << std::endl;
              for (unsigned int pt_idx = 0; pt_idx < output_points[c].size();
                   pt_idx++)
                deallog << output_points[c][pt_idx] << " from process "
                        << output_ranks[c][pt_idx] << " to " << my_rank
                        << std::endl;
              deallog << "Expected points:" << std::endl;
              for (unsigned int pt_idx = 0;
                   pt_idx < computed_points[c_cell].size();
                   pt_idx++)
                deallog << computed_points[c_cell][pt_idx] << std::endl;
            }
          else
            {
              // Checking if the points inside are the same
              for (unsigned int pt_idx = 0; pt_idx < output_points[c].size();
                   pt_idx++)
                {
                  const auto &pt    = output_points[c][pt_idx];
                  auto        pt_it = std::find(computed_points[c_cell].begin(),
                                         computed_points[c_cell].end(),
                                         pt);
                  if (pt_it == computed_points[c_cell].end())
                    {
                      deallog << "ERROR: point " << pt << " not found"
                              << std::endl;
                      test_passed = false;
                    }
                  else
                    {
                      unsigned int c_pt =
                        pt_it - computed_points[c_cell].begin();
                      // Checking the value of the transformed point
                      if ((output_qpoints[c][pt_idx] -
                           computed_qpoints[c_cell][c_pt])
                            .norm() > 1e-12)
                        {
                          // Cell not found: adding a new cell
                          deallog << "ERROR: qpoint " << c_pt << " not matching"
                                  << std::endl;
                          test_passed = false;
                        }
                      // Checking the rank of the owner
                      if (output_ranks[c][pt_idx] !=
                          computed_ranks[c_cell][c_pt])
                        {
                          // Cell not found: adding a new cell
                          deallog << "ERROR: rank of point " << c_pt
                                  << " not matching" << std::endl;
                          test_passed = false;
                        }
                    }
                }
            }
        }
    }



  if (output_computed_pts != computed_pts)
    {
      deallog << "ERROR: the number of points is different from expected: "
              << std::endl;
      deallog << "Number of locally computed points: " << computed_pts
              << std::endl;
      deallog << "Number of points from distributed: " << output_computed_pts
              << std::endl;
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

  deallog << "Deal.II distributed_compute_point_locations:" << std::endl;
  deallog << "Test on parallel setting 2D:" << std::endl;
  test_compute_pt_loc<2>(3, 3);
  deallog << "Test on parallel setting 3D:" << std::endl;
  test_compute_pt_loc<3>(3, 2);
}
