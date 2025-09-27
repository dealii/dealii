// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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

#include <deal.II/matrix_free/portable_fe_evaluation.h>
#include <deal.II/matrix_free/portable_matrix_free.h>
#include <deal.II/matrix_free/tools.h>
#include <deal.II/matrix_free/vector_access_internal.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

#include <deal.II/numerics/matrix_creator.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim, int fe_degree, int n_components, typename Number>
class LaplaceOperatorQuad
{
public:
  DEAL_II_HOST_DEVICE void
  operator()(
    Portable::FEEvaluation<dim, fe_degree, fe_degree + 1, n_components, Number>
             *fe_eval,
    const int q_point) const
  {
    fe_eval->submit_gradient(fe_eval->get_gradient(q_point), q_point);
  }

  static const unsigned int n_q_points =
    dealii::Utilities::pow(fe_degree + 1, dim);

  static const unsigned int n_local_dofs =
    dealii::Utilities::pow(fe_degree + 1, dim);
};


template <int dim,
          int fe_degree,
          int n_points,
          int n_components,
          typename Number>
class Test
{
  using VectorType =
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>;

public:
  Test(const Portable::MatrixFree<dim, Number> &matrix_free,
       const AffineConstraints<Number>         &constraints)
    : matrix_free(matrix_free)
    , constraints(constraints)
  {}

  void
  do_test()
  {
    // compute diagonal with the new function
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>
      diagonal_global;
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>
      diagonal_global_host;
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>
      diagonal_global_reference;

    SparseMatrix<Number> A_ref;
    SparsityPattern      sparsity_pattern;

    {
      matrix_free.initialize_dof_vector(diagonal_global);
      LaplaceOperatorQuad<dim, fe_degree, n_components, Number>
        laplace_operator_quad;
      MatrixFreeTools::
        compute_diagonal<dim, fe_degree, n_points, n_components, Number>(
          matrix_free,
          diagonal_global,
          laplace_operator_quad,
          EvaluationFlags::gradients,
          EvaluationFlags::gradients);

      matrix_free.initialize_dof_vector(diagonal_global_host);
      LinearAlgebra::ReadWriteVector<Number> rw_vector(
        diagonal_global.get_partitioner()->locally_owned_range());
      rw_vector.import_elements(diagonal_global, VectorOperation::insert);
      diagonal_global_host.import_elements(rw_vector, VectorOperation::insert);
      for (unsigned int i = 0; i < diagonal_global.size(); ++i)
        {
          if (diagonal_global.get_partitioner()->in_local_range(i))
            {
              Assert(diagonal_global_host[i] > 0,
                     ExcMessage("Diagonal non-positive at position " +
                                std::to_string(i)));
            }
        }

      diagonal_global_host.print(deallog.get_file_stream());
    }

    const bool test_matrix =
      (Utilities::MPI::job_supports_mpi() == false) ||
      (Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1);

    if (test_matrix)
      {
        DynamicSparsityPattern dsp(matrix_free.get_dof_handler().n_dofs());
        DoFTools::make_sparsity_pattern(matrix_free.get_dof_handler(),
                                        dsp,
                                        constraints);
        sparsity_pattern.copy_from(dsp);
        A_ref.reinit(sparsity_pattern);

        Function<dim, Number> *scaling = nullptr;
        MatrixCreator::create_laplace_matrix(matrix_free.get_dof_handler(),
                                             QGauss<dim>(fe_degree + 1),
                                             A_ref,
                                             scaling,
                                             constraints);

        for (unsigned int i = 0; i < diagonal_global.size(); ++i)
          {
            if (!constraints.is_constrained(i))
              {
                Assert(std::abs(A_ref(i, i) - diagonal_global_host(i)) < 1e-6,
                       ExcMessage("Wrong diagonal entry at position " +
                                  std::to_string(i)));
              }
          }
      }
  }

  const Portable::MatrixFree<dim, Number> &matrix_free;
  const AffineConstraints<Number>         &constraints;
};
