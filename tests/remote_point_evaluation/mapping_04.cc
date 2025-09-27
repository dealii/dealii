// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test Utilities::MPI::RemotePointEvaluation if points have been found
// and if the mapping is unique.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools_evaluate.h>

#include "../tests.h"


template <int dim>
void
do_test(const bool                     enforce_unique_map,
        const std::vector<Point<dim>> &evaluation_points,
        const std::pair<bool, bool>   &expected_result)
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::subdivided_hyper_rectangle(tria, {2, 1}, {0, 0}, {2, 1});

  MappingQ1<dim> mapping;

  typename Utilities::MPI::RemotePointEvaluation<dim>::AdditionalData
    additional_data;
  additional_data.enforce_unique_mapping = enforce_unique_map;

  Utilities::MPI::RemotePointEvaluation<dim> eval(additional_data);

  eval.reinit(evaluation_points, tria, mapping);

  Assert(eval.is_map_unique() == expected_result.first, ExcInternalError());
  Assert(eval.all_points_found() == expected_result.second, ExcInternalError());
}


template <int dim>
void
test()
{
  const auto my_rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  {
    // (-, -)
    std::vector<Point<dim>> evaluation_points;

    do_test(false, evaluation_points, {true, true});
  }

  {
    // (not found, -)
    std::vector<Point<dim>> evaluation_points;

    if (my_rank == 0)
      evaluation_points.emplace_back(2.5, 0.5);

    do_test(false, evaluation_points, {false, false});
  }

  {
    // (unique, -)
    std::vector<Point<dim>> evaluation_points;

    if (my_rank == 0)
      evaluation_points.emplace_back(1.5, 0.5);

    do_test(false, evaluation_points, {true, true});
  }

  {
    // (non-unique, -)
    std::vector<Point<dim>> evaluation_points;

    if (my_rank == 0)
      evaluation_points.emplace_back(1.0, 0.5);

    do_test(false, evaluation_points, {false, true});
  }

  {
    // (enforced-unique, -)
    std::vector<Point<dim>> evaluation_points;

    if (my_rank == 0)
      evaluation_points.emplace_back(1.0, 0.5);

    do_test(true, evaluation_points, {true, true});
  }

  {
    // (-, not found)
    std::vector<Point<dim>> evaluation_points;

    if (my_rank == 1)
      evaluation_points.emplace_back(2.5, 0.5);

    do_test(false, evaluation_points, {false, false});
  }

  {
    // (-, unique)
    std::vector<Point<dim>> evaluation_points;

    if (my_rank == 1)
      evaluation_points.emplace_back(0.5, 0.5);

    do_test(false, evaluation_points, {true, true});
  }

  {
    // (-, non-unique)
    std::vector<Point<dim>> evaluation_points;

    if (my_rank == 1)
      evaluation_points.emplace_back(1.0, 0.5);

    do_test(false, evaluation_points, {false, true});
  }

  {
    // (-, enforced-unique)
    std::vector<Point<dim>> evaluation_points;

    if (my_rank == 1)
      evaluation_points.emplace_back(1.0, 0.5);

    do_test(true, evaluation_points, {true, true});
  }

  {
    // (unique, not found)
    std::vector<Point<dim>> evaluation_points;

    if (my_rank == 0)
      evaluation_points.emplace_back(1.5, 0.5);
    if (my_rank == 1)
      evaluation_points.emplace_back(2.5, 0.5);

    do_test(false, evaluation_points, {false, false});
  }

  {
    // (unique, not unique)
    std::vector<Point<dim>> evaluation_points;

    if (my_rank == 0)
      evaluation_points.emplace_back(1.5, 0.5);
    if (my_rank == 1)
      evaluation_points.emplace_back(1.0, 0.5);

    do_test(false, evaluation_points, {false, true});
  }

  {
    // (unique, enforced unique)
    std::vector<Point<dim>> evaluation_points;

    if (my_rank == 0)
      evaluation_points.emplace_back(1.5, 0.5);
    if (my_rank == 1)
      evaluation_points.emplace_back(1.0, 0.5);

    do_test(true, evaluation_points, {true, true});
  }
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    all;

  AssertDimension(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD), 2);

  test<2>();
}
