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

// Test Utilities::MPI::RemotePointEvaluation with and without unique mapping.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools_evaluate.h>

#include "../tests.h"



template <int dim>
void
test(const bool enforce_unique_map)
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

  std::vector<Point<dim>> evaluation_points;

  for (unsigned int j = 0; j <= 4; ++j)
    for (unsigned int i = 0; i <= 4; ++i)
      evaluation_points.emplace_back(i * 0.25, j * 0.25);

  typename Utilities::MPI::RemotePointEvaluation<dim>::AdditionalData
    additional_data;
  additional_data.enforce_unique_mapping = enforce_unique_map;
  additional_data.tolerance              = 1e-6;
  Utilities::MPI::RemotePointEvaluation<dim> eval(additional_data);

  const auto result_avg =
    VectorTools::point_values<1>(mapping,
                                 dof_handler,
                                 vec,
                                 evaluation_points,
                                 eval,
                                 VectorTools::EvaluationFlags::avg);
  const auto result_min = VectorTools::point_values<1>(
    eval, dof_handler, vec, VectorTools::EvaluationFlags::min);
  const auto result_max = VectorTools::point_values<1>(
    eval, dof_handler, vec, VectorTools::EvaluationFlags::max);
  const auto result_insert = VectorTools::point_values<1>(
    eval, dof_handler, vec, VectorTools::EvaluationFlags::insert);

  for (unsigned int i = 0; i < evaluation_points.size(); ++i)
    deallog << i << ' ' << result_avg[i] << ' ' << result_min[i] << ' '
            << result_max[i] << ' ' << result_insert[i] << ' '
            << eval.get_point_ptrs()[i + 1] - eval.get_point_ptrs()[i]
            << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    all;

  test<2>(false);
  test<2>(true);
}
