// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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
          int n_points,
          typename Number,
          typename VectorizedArrayType>
class Tester
{
public:
  void
  cell_function(
    std::vector<
      FEEvaluation<dim, fe_degree, n_points, 1, Number, VectorizedArrayType>>
      &phi) const
  {
    this->cell_operation(phi);
  }

  std::function<void(
    std::vector<
      FEEvaluation<dim, fe_degree, n_points, 1, Number, VectorizedArrayType>>
      &)>
    cell_operation;
};

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
  // DoFTools::make_hanging_node_constraints(dof_handler_system,
  //                                        constraints_system);
  // DoFTools::make_zero_boundary_constraints(dof_handler_system,
  //                                         constraints_system);
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
    mapping, dof_handler_system, constraints_system, quad, additional_data);

  const auto kernel =
    [](std::vector<
       FEEvaluation<dim, fe_degree, n_points, 1, Number, VectorizedArrayType>>
         &phi) {
      for (unsigned int c = 0; c < phi.size(); ++c)
        phi[c].evaluate(true, false, false);

      for (unsigned int q = 0; q < phi[0].n_q_points; ++q)
        {
          std::vector<VectorizedArrayType> values(n_components);
          std::vector<VectorizedArrayType> values_temp(n_components);

          for (unsigned int c = 0; c < phi.size(); ++c)
            {
              values[c]      = phi[c].get_value(q);
              values_temp[c] = phi[c].get_value(q);
            }

          for (unsigned int c = 0; c < n_components; ++c)
            {
              values[c] *= (c + 1.0);

              for (unsigned int c1 = c + 1; c1 < n_components; ++c1)
                values[c] += values_temp[c1] * (c + 1.0);
            }

          for (unsigned int c = 0; c < phi.size(); ++c)
            phi[c].submit_value(values[c], q);
        }

      for (unsigned int c = 0; c < phi.size(); ++c)
        phi[c].integrate(true, false);
    };

  Tester<dim, fe_degree, n_points, Number, VectorizedArrayType> test;
  test.cell_operation = kernel;

  const auto kernel_vec = [](FEEvaluation<dim,
                                          fe_degree,
                                          n_points,
                                          n_components,
                                          Number,
                                          VectorizedArrayType> &phi) {
    phi.evaluate(true, false, false);
    for (unsigned int q = 0; q < phi.n_q_points; ++q)
      {
        auto values      = phi.get_value(q);
        auto values_temp = phi.get_value(q);

        for (unsigned int c = 0; c < n_components; ++c)
          {
            values[c] *= (c + 1.0);

            for (unsigned int c1 = c + 1; c1 < n_components; ++c1)
              values[c] += values_temp[c1] * (c + 1.0);
          }

        phi.submit_value(values, q);
      }
    phi.integrate(true, false);
  };

  {
    LinearAlgebra::distributed::Vector<Number> diagonal_global_1;
    matrix_free.initialize_dof_vector(diagonal_global_1);

    MatrixFreeTools::compute_diagonal<dim,
                                      fe_degree,
                                      n_points,
                                      n_components,
                                      Number,
                                      VectorizedArrayType>(matrix_free,
                                                           diagonal_global_1,
                                                           kernel_vec);

    diagonal_global_1.print(deallog.get_file_stream());

    LinearAlgebra::distributed::Vector<Number> diagonal_global_2;
    matrix_free.initialize_dof_vector(diagonal_global_2);

    MatrixFreeTools::
      compute_diagonal<dim, fe_degree, n_points, Number, VectorizedArrayType>(
        matrix_free, diagonal_global_2, n_components, kernel);

    diagonal_global_2.print(deallog.get_file_stream());

    LinearAlgebra::distributed::Vector<Number> diagonal_global_3;
    matrix_free.initialize_dof_vector(diagonal_global_3);

    MatrixFreeTools::compute_diagonal(
      matrix_free,
      diagonal_global_3,
      n_components,
      &Tester<dim, fe_degree, n_points, Number, VectorizedArrayType>::
        cell_function,
      &test);

    diagonal_global_3.print(deallog.get_file_stream());
  }


  {
    SparseMatrix<Number> A1, A2, A3;
    SparsityPattern      sparsity_pattern;

    DynamicSparsityPattern dsp(matrix_free.get_dof_handler().n_dofs());
    DoFTools::make_sparsity_pattern(matrix_free.get_dof_handler(),
                                    dsp,
                                    constraints_system);
    sparsity_pattern.copy_from(dsp);
    A1.reinit(sparsity_pattern);
    A2.reinit(sparsity_pattern);
    A3.reinit(sparsity_pattern);

    MatrixFreeTools::compute_matrix<dim,
                                    fe_degree,
                                    n_points,
                                    Number,
                                    VectorizedArrayType,
                                    SparseMatrix<Number>>(
      matrix_free, constraints_system, A1, n_components, kernel);

    MatrixFreeTools::compute_matrix<dim,
                                    fe_degree,
                                    n_points,
                                    n_components,
                                    Number,
                                    VectorizedArrayType,
                                    SparseMatrix<Number>>(matrix_free,
                                                          constraints_system,
                                                          A2,
                                                          kernel_vec);

    MatrixFreeTools::compute_matrix(
      matrix_free,
      constraints_system,
      A3,
      n_components,
      &Tester<dim, fe_degree, n_points, Number, VectorizedArrayType>::
        cell_function,
      &test);

    for (auto entry_1 = A1.begin(), entry_2 = A2.begin(), entry_3 = A3.begin();
         entry_1 != A1.end();
         ++entry_1, ++entry_2, ++entry_3)
      {
        deallog << entry_1->row() << " " << entry_1->column() << " "
                << entry_1->value() << " " << entry_2->value() << " "
                << entry_3->value() << std::endl;
        Assert(std::abs(entry_1->value() - entry_2->value()) < 1e-5,
               ExcInternalError());
        Assert(std::abs(entry_1->value() - entry_3->value()) < 1e-5,
               ExcInternalError());
      }

    deallog << "OK!" << std::endl;
  }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test<2, 1>();
}
