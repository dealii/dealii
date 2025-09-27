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

// Test GridTools::distributed_compute_point_locations for the serial case

#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
test_compute_pt_loc(unsigned int n_points)
{
  MPI_Comm mpi_communicator = MPI_COMM_WORLD;
  deallog << "Testing for dim = " << dim << std::endl;
  deallog << "Testing on: " << n_points << " points." << std::endl;

  // Creating a grid in the square [0,1]x[0,1]
  parallel::distributed::Triangulation<dim> tria(mpi_communicator);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(std::max(6 - dim, 2));

  // Creating the random points
  std::vector<Point<dim>> points;

  for (std::size_t i = 0; i < n_points; ++i)
    points.push_back(random_point<dim>());

  // Initializing the cache
  GridTools::Cache<dim, dim> cache(tria);

  // Computing the description of the locally owned part of the mesh
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
    distributed_compute_point_locations(cache, points, global_bboxes);
  // Testing in serial against the serial version
  auto cell_qpoint_map = GridTools::compute_point_locations(cache, points);

  const auto &serial_cells   = std::get<0>(cell_qpoint_map);
  const auto &serial_qpoints = std::get<1>(cell_qpoint_map);
  const auto &serial_map     = std::get<2>(cell_qpoint_map);
  std::size_t n_cells        = std::get<0>(output_tuple).size();

  deallog << "Points found in " << n_cells << " cells" << std::endl;

  // testing if the result coincides with
  // the serial one
  for (unsigned int c = 0; c < n_cells; ++c)
    {
      const auto &cell            = std::get<0>(output_tuple)[c];
      const auto &quad            = std::get<1>(output_tuple)[c];
      const auto &local_map       = std::get<2>(output_tuple)[c];
      const auto &original_points = std::get<3>(output_tuple)[c];
      const auto &ranks           = std::get<4>(output_tuple)[c];

      const auto pos_cell =
        std::find(serial_cells.begin(), serial_cells.end(), cell);
      for (auto r : ranks)
        if (r != 0)
          deallog << "ERROR: rank is not 0 but " << std::to_string(r)
                  << std::endl;

      if (pos_cell == serial_cells.end())
        deallog << "ERROR: cell not found" << std::endl;
      else
        {
          auto serial_cell_idx = pos_cell - serial_cells.begin();
          if (original_points.size() != serial_qpoints[serial_cell_idx].size())
            deallog << "ERROR: in the number of points for cell"
                    << std::to_string(serial_cell_idx) << std::endl;
          if (quad.size() != serial_qpoints[serial_cell_idx].size())
            deallog << "ERROR: in the number of points for cell"
                    << std::to_string(serial_cell_idx) << std::endl;

          unsigned int pt_num           = 0;
          const auto   serial_local_map = serial_map[serial_cell_idx];
          for (const auto &p_idx : local_map)
            {
              auto serial_pt_pos = std::find(serial_local_map.begin(),
                                             serial_local_map.end(),
                                             p_idx);
              auto serial_pt_idx = serial_pt_pos - serial_local_map.begin();
              if (serial_pt_pos == serial_local_map.end())
                deallog << "ERROR: point index not found for "
                        << std::to_string(serial_pt_idx) << std::endl;
              else
                {
                  if ((original_points[pt_num] - points[p_idx]).norm() > 1e-12)
                    {
                      deallog
                        << "ERROR: Point in serial : " << points[p_idx]
                        << " Point in distributed: " << original_points[pt_num]
                        << std::endl;
                    }

                  if ((quad[pt_num] -
                       serial_qpoints[serial_cell_idx][serial_pt_idx])
                        .norm() > 1e-10)
                    {
                      deallog
                        << " ERROR: Transformation of qpoint to point is not correct"
                        << std::endl;
                      deallog << "qpoint in serial : " << quad[pt_num]
                              << " Point in distributed: "
                              << serial_qpoints[serial_cell_idx][serial_pt_idx]
                              << std::endl;
                    }
                }
              ++pt_num;
            }
        }
    }

  deallog << "Test finished" << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog << "Deal.II distributed_compute_point_locations:" << std::endl;
  test_compute_pt_loc<2>(100);
  test_compute_pt_loc<3>(200);
}
