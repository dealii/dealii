// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test Portable::MatrixFree with level cells. Creates a hypercube grid,
// refines it 3 times, distributes dofs and mg dofs. Then initializes
// Portable MatrixFree and uses it to evaluate Laplace operator applied
// to a random vector on level 0. Compares with classical matrix-based
// evaluation.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/portable_fe_evaluation.h>
#include <deal.II/matrix_free/portable_matrix_free.h>

#include <deal.II/multigrid/mg_tools.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

#include "Kokkos_Core.hpp"

using namespace dealii;


template <int dim, int fe_degree, typename Number, int n_q_points_1d>
class LaplaceOperatorQuad
{
public:
  DEAL_II_HOST_DEVICE
  LaplaceOperatorQuad()
  {}

  DEAL_II_HOST_DEVICE void
  operator()(
    Portable::FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> *phi,
    const int                                                         q) const
  {
    phi->submit_gradient(phi->get_gradient(q), q);
  }

  static const unsigned int n_q_points =
    dealii::Utilities::pow(n_q_points_1d, dim);
};



template <int dim, int fe_degree, typename Number, int n_q_points_1d>
class LaplaceOperator
{
public:
  static const unsigned int n_q_points =
    dealii::Utilities::pow(n_q_points_1d, dim);

  LaplaceOperator()
  {}

  DEAL_II_HOST_DEVICE void
  operator()(const typename Portable::MatrixFree<dim, Number>::Data *data,
             const Portable::DeviceVector<Number>                   &src,
             Portable::DeviceVector<Number>                         &dst) const
  {
    Portable::FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> fe_eval(
      data);
    fe_eval.read_dof_values(src);
    fe_eval.evaluate(EvaluationFlags::gradients);

    LaplaceOperatorQuad<dim, fe_degree, Number, n_q_points_1d> quad;
    data->for_each_quad_point(
      [&](const int &q_point) { quad(&fe_eval, q_point); });

    fe_eval.integrate(EvaluationFlags::gradients);
    fe_eval.distribute_local_to_global(dst);
  }
};



template <int dim, int fe_degree, typename Number>
void
test()
{
  const unsigned int test_level    = 1;
  const unsigned int n_q_points_1d = fe_degree + 1;

  // Create triangulation with mesh smoothing for multigrid
  Triangulation<dim> tria(
    Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);

  // Setup DoFHandler
  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();

  // Constraints
  AffineConstraints<Number> constraints;
  constraints.close();

  deallog << "Testing Portable::MatrixFree on level " << test_level
          << std::endl;
  deallog << "Number of dofs on level " << test_level << ": "
          << dof_handler.n_dofs(test_level) << std::endl;

  // Initialize Portable::MatrixFree for level cells
  MappingQ<dim>                                              mapping(fe_degree);
  Portable::MatrixFree<dim, Number>                          mf_data;
  typename Portable::MatrixFree<dim, Number>::AdditionalData additional_data;
  additional_data.mapping_update_flags = update_gradients | update_JxW_values;
  additional_data.mg_level             = test_level;
  const QGauss<1> quad(n_q_points_1d);

  mf_data.reinit(mapping, dof_handler, constraints, quad, additional_data);

  // Check that mg_level is stored correctly
  AssertThrow(mf_data.get_mg_level() == test_level,
              ExcMessage("mg_level not correctly stored"));

  // Create test vectors
  const unsigned int n_dofs = dof_handler.n_dofs(test_level);

  LinearAlgebra::distributed::Vector<Number, MemorySpace::Default> src(n_dofs);
  LinearAlgebra::distributed::Vector<Number, MemorySpace::Default> dst_device(
    n_dofs);
  Vector<Number> src_host(n_dofs);
  Vector<Number> dst_host(n_dofs);

  // Fill source vector with random values
  LinearAlgebra::ReadWriteVector<Number> src_rw(n_dofs);
  for (unsigned int i = 0; i < n_dofs; ++i)
    {
      if (constraints.is_constrained(i))
        continue;
      const double entry = Testing::rand() / (double)RAND_MAX;
      src_rw(i)          = entry;
      src_host(i)        = entry;
    }

  src.import_elements(src_rw, VectorOperation::insert);

  // Apply operator using Portable::MatrixFree
  mf_data.cell_loop(LaplaceOperator<dim, fe_degree, Number, n_q_points_1d>(),
                    src,
                    dst_device);
  Kokkos::fence();

  // Copy result back to host
  LinearAlgebra::ReadWriteVector<Number> dst_rw(n_dofs);
  dst_rw.import_elements(dst_device, VectorOperation::insert);

  // Assemble sparse matrix for comparison on the given level
  SparsityPattern sparsity;
  {
    DynamicSparsityPattern dsp(n_dofs, n_dofs);
    MGTools::make_sparsity_pattern(dof_handler, dsp, test_level);
    sparsity.copy_from(dsp);
  }

  SparseMatrix<double> sparse_matrix(sparsity);

  // Manually assemble the Laplace matrix on the level using FEValues
  {
    QGauss<dim>   quadrature_formula(fe_degree + 1);
    FEValues<dim> fe_values(mapping,
                            fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.cell_iterators_on_level(test_level))
      {
        cell_matrix = 0;
        fe_values.reinit(cell);

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          {
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  cell_matrix(i, j) +=
                    (fe_values.shape_grad(i, q_point) *
                     fe_values.shape_grad(j, q_point) * fe_values.JxW(q_point));
              }
          }

        cell->get_mg_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix,
                                               local_dof_indices,
                                               sparse_matrix);
      }
  }

  // Apply sparse matrix
  sparse_matrix.vmult(dst_host, src_host);

  // Compare results
  Number error_norm = 0.;
  for (unsigned int i = 0; i < n_dofs; ++i)
    error_norm += std::pow(dst_rw(i) - dst_host(i), 2);

  const double diff_norm = std::sqrt(error_norm) / dst_host.linfty_norm();

  if (diff_norm > 1e-12)
    {
      deallog << "Error: Norm of difference " << diff_norm
              << " exceeds tolerance!" << std::endl;
      AssertThrow(false, ExcMessage("Result does not match reference!"));
    }
  else
    {
      deallog << "Norm of difference OK " << std::endl << std::endl;
    }
}


int
main()
{
  initlog();
  deallog.depth_console(0);

  deallog << std::setprecision(8);

  Kokkos::initialize();

  {
    deallog.push("2d");
    test<2, 1, double>();
    test<2, 2, double>();
    deallog.pop();
  }

  Kokkos::finalize();

  return 0;
}
