// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/tools.h>
#include <deal.II/matrix_free/vector_access_internal.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

using namespace dealii;


template <int dim,
          int fe_degree,
          int n_points,
          int n_components,
          typename Number,
          typename VectorizedArrayType>
class Test
{
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

public:
  Test(const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
       const AffineConstraints<Number> &                   constraints,
       const std::function<void(FEEvaluation<dim,
                                             fe_degree,
                                             n_points,
                                             n_components,
                                             Number,
                                             VectorizedArrayType> &)>
         &cell_operation)
    : matrix_free(matrix_free)
    , constraints(constraints)
    , cell_operation(cell_operation)
  {}

  void
  do_test()
  {
    // compute diagonal with the new function
    VectorType diagonal_global, diagonal_global_reference;

    SparseMatrix<Number> A1, A2, A_ref;
    SparsityPattern      sparsity_pattern;

    const bool test_matrix = Utilities::MPI::n_mpi_processes(
                               matrix_free.get_task_info().communicator) == 1;

    if (test_matrix)
      {
        DynamicSparsityPattern dsp(matrix_free.get_dof_handler().n_dofs());
        DoFTools::make_sparsity_pattern(matrix_free.get_dof_handler(),
                                        dsp,
                                        constraints);
        sparsity_pattern.copy_from(dsp);
        A1.reinit(sparsity_pattern);
        A2.reinit(sparsity_pattern);
        A_ref.reinit(sparsity_pattern);
      }

    {
      MatrixFreeTools::compute_diagonal<dim,
                                        fe_degree,
                                        n_points,
                                        n_components,
                                        Number,
                                        VectorizedArrayType>(
        matrix_free, diagonal_global, [&](auto &phi) {
          this->cell_operation(phi);
        });

      diagonal_global.print(deallog.get_file_stream());
      deallog << diagonal_global.l2_norm() << std::endl;
    }

    {
      VectorType diagonal_global;
      MatrixFreeTools::compute_diagonal(matrix_free,
                                        diagonal_global,
                                        &Test::cell_function,
                                        this);

      diagonal_global.print(deallog.get_file_stream());
      deallog << diagonal_global.l2_norm() << std::endl;
    }

    if (test_matrix)
      {
        MatrixFreeTools::compute_matrix<dim,
                                        fe_degree,
                                        n_points,
                                        n_components,
                                        Number,
                                        VectorizedArrayType,
                                        SparseMatrix<Number>>(
          matrix_free, constraints, A1, [&](auto &phi) {
            this->cell_operation(phi);
          });
      }

    if (test_matrix)
      {
        MatrixFreeTools::compute_matrix(
          matrix_free, constraints, A2, &Test::cell_function, this);
      }

    // compute diagonal globally
    {
      VectorType src, temp;

      matrix_free.initialize_dof_vector(src);
      matrix_free.initialize_dof_vector(diagonal_global_reference);
      matrix_free.initialize_dof_vector(temp);

      for (unsigned int i = 0; i < src.size(); i++)
        {
          if (src.get_partitioner()->in_local_range(i))
            src[i] = 1.0;

          matrix_free.cell_loop(&Test::cell_operation_range, this, temp, src);

          if (src.get_partitioner()->in_local_range(i))
            {
              diagonal_global_reference[i] = temp[i];
              src[i]                       = 0.0;
            }

          if (test_matrix)
            {
              for (unsigned int j = 0; j < src.size(); j++)
                if (temp[j] != 0.0)
                  A_ref(j, i) = temp[j];
                else if (i == j)
                  A_ref(j, i) = 1.0;
            }

          temp = 0.0;
        }

      diagonal_global_reference.print(deallog.get_file_stream());
      deallog << diagonal_global_reference.l2_norm() << std::endl;
    }

    if (test_matrix)
      {
        auto a1    = A1.begin();
        auto a2    = A2.begin();
        auto a_ref = A_ref.begin();

        for (; a1 != A1.end(); ++a1, ++a2, ++a_ref)
          {
            if (a1->row() == a1->column() &&
                constraints.is_constrained(a1->row()))
              continue;

            Assert(std::abs(a1->value() - a_ref->value()) < 1e-6,
                   ExcNotImplemented());
            Assert(std::abs(a2->value() - a_ref->value()) < 1e-6,
                   ExcNotImplemented());
          }
      }
  }

  void
  cell_operation_range(const MatrixFree<dim, Number, VectorizedArrayType> &data,
                       VectorType &                                        dst,
                       const VectorType &                                  src,
                       const std::pair<unsigned int, unsigned int> &pair) const
  {
    FEEvaluation<dim,
                 fe_degree,
                 n_points,
                 n_components,
                 Number,
                 VectorizedArrayType>
      phi(data, pair);
    for (auto cell = pair.first; cell < pair.second; cell++)
      {
        phi.reinit(cell);
        phi.read_dof_values(src);
        this->cell_operation(phi);
        phi.distribute_local_to_global(dst);
      }
  }

  void
  cell_function(FEEvaluation<dim,
                             fe_degree,
                             n_points,
                             n_components,
                             Number,
                             VectorizedArrayType> &phi) const
  {
    this->cell_operation(phi);
  }

  const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free;
  const AffineConstraints<Number> &                   constraints;
  const std::function<void(FEEvaluation<dim,
                                        fe_degree,
                                        n_points,
                                        n_components,
                                        Number,
                                        VectorizedArrayType> &)>
    cell_operation;
};
