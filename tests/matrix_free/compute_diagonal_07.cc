// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2022 by the deal.II authors
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



// Similar to compute_diagonal_02 but testing block vectorss.

#include <deal.II/lac/la_parallel_block_vector.h>

#include "compute_diagonal_util.h"

using namespace dealii;

template <int dim,
          int fe_degree,
          int n_points                 = fe_degree + 1,
          int n_components             = dim,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_ball(tria);
  tria.refine_global(0);

  const FE_Q<dim>     fe_q(fe_degree);
  const FESystem<dim> fe_system(fe_q, n_components);

  DoFHandler<dim> dof_handler_q(tria);
  dof_handler_q.distribute_dofs(fe_q);

  DoFHandler<dim> dof_handler_system(tria);
  dof_handler_system.distribute_dofs(fe_system);

  AffineConstraints<Number> constraints_q;
  DoFTools::make_hanging_node_constraints(dof_handler_q, constraints_q);
  DoFTools::make_zero_boundary_constraints(dof_handler_q, constraints_q);
  constraints_q.close();

  AffineConstraints<Number> constraints_system;
  DoFTools::make_hanging_node_constraints(dof_handler_system,
                                          constraints_system);
  DoFTools::make_zero_boundary_constraints(dof_handler_system,
                                           constraints_system);
  constraints_system.close();

  typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
    additional_data;
  additional_data.mapping_update_flags = update_values | update_gradients;
  additional_data.tasks_parallel_scheme =
    MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData::
      TasksParallelScheme::none;

  MappingQ<dim> mapping(1);
  QGauss<1>     quad(fe_degree + 1);

  MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
  matrix_free.reinit(
    mapping,
    std::vector<const DoFHandler<dim> *>{&dof_handler_system, &dof_handler_q},
    std::vector<const AffineConstraints<Number> *>{&constraints_system,
                                                   &constraints_q},
    quad,
    additional_data);

  const auto kernel = [](FEEvaluation<dim,
                                      fe_degree,
                                      n_points,
                                      n_components,
                                      Number,
                                      VectorizedArrayType> &phi) {
    phi.evaluate(false, true, false);
    for (unsigned int q = 0; q < phi.n_q_points; ++q)
      {
        phi.submit_symmetric_gradient(2.0 * phi.get_symmetric_gradient(q), q);
      }
    phi.integrate(false, true);
  };

  LinearAlgebra::distributed::Vector<Number> diagonal_global;
  matrix_free.initialize_dof_vector(diagonal_global, 0);

  MatrixFreeTools::compute_diagonal<dim,
                                    fe_degree,
                                    n_points,
                                    n_components,
                                    Number,
                                    VectorizedArrayType>(matrix_free,
                                                         diagonal_global,
                                                         kernel,
                                                         0);

  diagonal_global.print(deallog.get_file_stream());

  LinearAlgebra::distributed::BlockVector<Number> diagonal_global_block(dim);

  for (unsigned int comp = 0; comp < dim; ++comp)
    matrix_free.initialize_dof_vector(diagonal_global_block.block(comp), 1);

  MatrixFreeTools::compute_diagonal<dim,
                                    fe_degree,
                                    n_points,
                                    n_components,
                                    Number,
                                    VectorizedArrayType>(matrix_free,
                                                         diagonal_global_block,
                                                         kernel,
                                                         1);

  diagonal_global_block.print(deallog.get_file_stream());
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test<2, 1>();
}
