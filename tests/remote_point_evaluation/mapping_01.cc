// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2022 by the deal.II authors
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

// Test Utilities::MPI::RemotePointEvaluation with and without unique mapping.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_vector.h>

#include <deal.II/numerics/vector_tools_evaluate.h>

#include "../tests.h"

using namespace dealii;


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

  Utilities::MPI::RemotePointEvaluation<dim> eval(1e-6, enforce_unique_map);

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
