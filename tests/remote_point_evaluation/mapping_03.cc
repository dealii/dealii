// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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

// Test Utilities::MPI::RemotePointEvaluation with and without unique mapping
// for a within a tolerance.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

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

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  Vector<double> vec(dof_handler.n_dofs());

  vec = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) + 1;

  MappingQ1<dim> mapping;

  std::vector<Point<dim>> evaluation_points;

  for (int i = 1; i < 15; ++i)
    // should be assigned to rank 0 for unique mapping
    evaluation_points.emplace_back(0.5, 0.5 + std::pow<double>(10.0, -i));

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
    deallog << "1e-" << (i + 1) << ' ' << result_avg[i] << ' ' << result_min[i]
            << ' ' << result_max[i] << ' ' << result_insert[i] << ' '
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
