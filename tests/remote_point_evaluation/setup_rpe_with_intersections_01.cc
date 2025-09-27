// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Face intersections 1D-2D
// Test RPE setup with distributed_compute_intersection_locations()

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_remote_point_evaluation.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



template <int dim>
std::vector<std::vector<BoundingBox<dim>>>
get_global_bboxes(const Triangulation<dim> &tria,
                  const Mapping<dim>       &mapping,
                  const unsigned int        rtree_level = 0)
{
  std::vector<dealii::BoundingBox<dim>> local_boxes;
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      local_boxes.emplace_back(mapping.get_bounding_box(cell));

  // create r-tree of bounding boxes
  const auto local_tree = pack_rtree(local_boxes);

  // compress r-tree to a minimal set of bounding boxes
  std::vector<std::vector<BoundingBox<dim>>> global_bboxes(1);
  global_bboxes[0] = extract_rtree_level(local_tree, rtree_level);

  return global_bboxes;
}

template <int structdim, int dim>
void
do_test(const unsigned int n_quad_points)
{
  // create triangulation
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  // create cache
  MappingQ1<dim>              mapping;
  const GridTools::Cache<dim> cache(tria, mapping);

  // create intersection requests
  std::vector<std::vector<Point<dim>>> intersection_requests;
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      Triangulation<dim> tria_other;
      // in this test we have to ensure quadrature points on intersections can
      // not found in different cells. This would lead to quadrature points
      // found twice in the intersections. Since we fetch the quadrature
      // points from found intersections an additional (not expected)
      // duplication takes place in RPE setup with these points.
      GridGenerator::hyper_cube(tria_other, 0.24, 0.76);
      const auto cell = tria_other.begin_active();

      for (const auto f : cell->face_indices())
        {
          const unsigned int      n_vertices = cell->face(f)->n_vertices();
          std::vector<Point<dim>> vertices(n_vertices);
          std::copy_n(mapping.get_vertices(cell, f).begin(),
                      n_vertices,
                      vertices.begin());

          intersection_requests.emplace_back(vertices);
        }
    }

  auto intersection_location =
    GridTools::internal::distributed_compute_intersection_locations<structdim>(
      cache,
      intersection_requests,
      get_global_bboxes<dim>(tria, mapping),
      std::vector<bool>(),
      1.0e-9);

  auto rpe_intersection_data =
    intersection_location
      .convert_to_distributed_compute_point_locations_internal(n_quad_points,
                                                               tria,
                                                               mapping);

  // setup rpe without additional search
  dealii::Utilities::MPI::RemotePointEvaluation<dim> rpe_intersections;
  rpe_intersections.reinit(rpe_intersection_data, tria, mapping);


  // get quadrature points of intersections
  std::vector<Point<dim>> points;

  const auto &recv_comps = intersection_location.recv_components;
  const QGaussSimplex<structdim> quadrature(n_quad_points);
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      for (const auto &rc : recv_comps)
        {
          const Quadrature<dim> &quad =
            quadrature.compute_affine_transformation(std::get<2>(rc));

          for (unsigned int q = 0; q < quad.size(); ++q)
            points.push_back(quad.point(q));
        }
    }

  // setup rpe using addition point search
  dealii::Utilities::MPI::RemotePointEvaluation<dim> rpe_points;
  rpe_points.reinit(points, tria, mapping);

  // compate rpe_points and rpe_intersections
  const auto &cell_data_inter = rpe_intersections.get_cell_data();
  const auto &cell_data_point = rpe_points.get_cell_data();

  // check CellData
  AssertThrow(cell_data_inter.cells.size() == cell_data_point.cells.size(),
              ExcMessage("CellData().cells.size() not OK"));
  deallog << "CellData().cells.size() OK" << std::endl;
  for (unsigned int i = 0; i < cell_data_inter.cells.size(); ++i)
    {
      AssertThrow(cell_data_inter.cells[i].first ==
                    cell_data_point.cells[i].first,
                  ExcMessage("CellData().cells.first not OK"));
      AssertThrow(cell_data_inter.cells[i].second ==
                    cell_data_point.cells[i].second,
                  ExcMessage("CellData().cells.second not OK"));
    }
  deallog << "CellData().cells OK" << std::endl;

  AssertThrow(cell_data_inter.reference_point_ptrs.size() ==
                cell_data_point.reference_point_ptrs.size(),
              ExcMessage("CellData().reference_point_ptrs.size() not OK"));
  deallog << "CellData().reference_point_ptrs.size() OK" << std::endl;
  for (unsigned int i = 0; i < cell_data_inter.reference_point_ptrs.size(); ++i)
    {
      AssertThrow(cell_data_inter.reference_point_ptrs[i] ==
                    cell_data_point.reference_point_ptrs[i],
                  ExcMessage("CellData().reference_point_ptrs not OK"));
    }
  deallog << "CellData().reference_point_ptrs OK" << std::endl;

  AssertThrow(cell_data_inter.reference_point_values.size() ==
                cell_data_point.reference_point_values.size(),
              ExcMessage("CellData().reference_point_values.size() not OK"));
  deallog << "CellData().reference_point_values.size() OK" << std::endl;
  for (unsigned int i = 0; i < cell_data_inter.reference_point_values.size();
       ++i)
    {
      AssertThrow(std::abs((cell_data_inter.reference_point_values[i] -
                            cell_data_point.reference_point_values[i])
                             .norm()) < 1.0e-9,
                  ExcMessage("CellData().reference_point_values not OK"));
    }
  deallog << "CellData().reference_point_values OK" << std::endl;
  deallog << "CellData() OK" << std::endl;

  // check get_point_ptrs()
  AssertThrow(rpe_intersections.get_point_ptrs().size() ==
                rpe_points.get_point_ptrs().size(),
              ExcMessage("get_point_ptrs().size() not OK"));
  deallog << "get_point_ptrs().size() OK" << std::endl;
  for (unsigned int i = 0; i < rpe_intersections.get_point_ptrs().size(); ++i)
    {
      AssertThrow(rpe_intersections.get_point_ptrs()[i] ==
                    rpe_points.get_point_ptrs()[i],
                  ExcMessage("get_point_ptrs() not OK"));
    }
  deallog << "get_point_ptrs() OK" << std::endl;

  // check is_map_unique()
  AssertThrow(rpe_intersections.is_map_unique() == rpe_points.is_map_unique(),
              ExcMessage("is_map_unique() not OK"));
  deallog << "is_map_unique() OK" << std::endl;

  // check all_points_found()
  AssertThrow(rpe_intersections.all_points_found() ==
                rpe_points.all_points_found(),
              ExcMessage("all_points_found() not OK"));
  deallog << "all_points_found() OK" << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;
  deallog.precision(8);

  for (unsigned int q = 1; q < 3; ++q)
    {
      deallog << "Face intersections 1D-2D with n_quadrature_points = " << q
              << std::endl;
      do_test<1, 2>(q);
    }
}
