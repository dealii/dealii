// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test GridTools::distributed_compute_point_locations() that it works with
// unique mapping.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/numerics/vector_tools_evaluate.h>

#include "../tests.h"



template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::subdivided_hyper_cube(tria, 2);
  tria.refine_global(1);

  FE_DGQ<dim>     fe(0);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  Vector<double> vec(dof_handler.n_dofs());

  for (const auto &cell : dof_handler.active_cell_iterators() |
                            IteratorFilters::LocallyOwnedCell())
    vec[cell->global_active_cell_index()] = cell->global_active_cell_index();

  MappingQ1<dim> mapping;

  std::vector<BoundingBox<dim>> local_boxes;
  for (const auto &cell : tria.cell_iterators_on_level(0))
    local_boxes.push_back(mapping.get_bounding_box(cell));

  // create r-tree of bounding boxes
  const auto local_tree = pack_rtree(local_boxes);

  // compress r-tree to a minimal set of bounding boxes
  const auto local_reduced_box = extract_rtree_level(local_tree, 0);

  // gather bounding boxes of other processes
  const auto global_bboxes =
    Utilities::MPI::all_gather(tria.get_mpi_communicator(), local_reduced_box);

  const GridTools::Cache<dim> cache(tria, mapping);

  std::vector<Point<dim>> local_points;

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    for (unsigned int j = 0; j <= 4; ++j)
      for (unsigned int i = 0; i <= 4; ++i)
        local_points.emplace_back(i * 0.25, j * 0.25);

  const auto temp =
    GridTools::distributed_compute_point_locations(cache,
                                                   local_points,
                                                   global_bboxes);

  unsigned int count = 0;

  for (const auto &i : std::get<1>(temp))
    count += i.size();

  const auto global_count = Utilities::MPI::sum(count, MPI_COMM_WORLD);
  const auto expected_global_count =
    Utilities::MPI::sum(local_points.size(), MPI_COMM_WORLD);

  deallog << global_count << ' ' << expected_global_count << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  initlog();

  test<2>();
}
