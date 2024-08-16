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

// Test VectorTools::point_values() for cell-data vectors.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_evaluate.h>

#include "../tests.h"


template <int dim, int spacedim>
std::shared_ptr<const Utilities::MPI::Partitioner>
create_partitioner(const DoFHandler<dim, spacedim> &dof_handler)
{
  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  return std::make_shared<const Utilities::MPI::Partitioner>(
    dof_handler.locally_owned_dofs(),
    locally_relevant_dofs,
    dof_handler.get_mpi_communicator());
}

void
test()
{
  const unsigned int dim             = 2;
  const unsigned int n_refinements_1 = 3;

  parallel::distributed::Triangulation<dim> tria_1(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria_1, -1.0, +1.0);
  tria_1.refine_global(n_refinements_1);
  DoFHandler<dim> dof_handler_1(tria_1);

  FE_DGQ<dim> fe_1(0);
  dof_handler_1.distribute_dofs(fe_1);

  // approach 1:
  LinearAlgebra::distributed::Vector<double> vector_1(
    create_partitioner(dof_handler_1));

  // approach 2:
  LinearAlgebra::distributed::Vector<double> vector_2(
    tria_1.global_active_cell_index_partitioner().lock());

  // approach 3:
  Vector<double> vector_3(tria_1.n_active_cells());

  for (const auto &cell : dof_handler_1.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        const auto value = cell->center()[0];

        // approach 1:
        Vector<double> temp(cell->get_fe().n_dofs_per_cell());
        temp = value;
        cell->set_dof_values(temp, vector_1);

        // approach 2:
        vector_2[cell->global_active_cell_index()] = value;

        // approach 3:
        vector_3[cell->active_cell_index()] = value;
      }

  const MappingQ1<dim> mapping_1;

  std::vector<Point<dim>> evaluation_points;

  const unsigned int n_intervals = 100;

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    for (unsigned int i = 0; i < n_intervals; ++i)
      {
        const double rad = 2.0 * numbers::PI / n_intervals * i;
        evaluation_points.emplace_back(0.75 * std::cos(rad),
                                       0.75 * std::sin(rad));
      }

  Utilities::MPI::RemotePointEvaluation<dim> evaluation_cache;
  const auto evaluation_point_results_1 = VectorTools::point_values<1>(
    mapping_1, dof_handler_1, vector_1, evaluation_points, evaluation_cache);
  const auto evaluation_point_results_2 =
    VectorTools::point_values<1>(evaluation_cache, tria_1, vector_2);
  const auto evaluation_point_results_3 =
    VectorTools::point_values<1>(evaluation_cache, tria_1, vector_3);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    for (unsigned int i = 0; i < n_intervals; ++i)
      {
        if (std::abs(evaluation_point_results_1[i] -
                     evaluation_point_results_2[i]) > 1e-10 ||
            std::abs(evaluation_point_results_1[i] -
                     evaluation_point_results_3[i]) > 1e-10)
          DEAL_II_NOT_IMPLEMENTED();
      }

  deallog << "OK!" << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  initlog();

  test();
}
