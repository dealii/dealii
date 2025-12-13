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



// Test Portable::MatrixFree with level cells and constraints.
// Applies Laplace operator to a random vector on the finest mg level and
// compares with result with using active cells.

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
  const unsigned int n_q_points_1d = fe_degree + 1;

  // Create triangulation with mesh smoothing for multigrid
  Triangulation<dim> tria(
    Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);

  const unsigned int test_level = tria.n_levels() - 1;

  // Setup DoFHandler
  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();

  MappingQ<dim> mapping(fe_degree);

  // Constraints
  AffineConstraints<Number> constraints;
  VectorTools::interpolate_boundary_values(
    mapping, dof_handler, 0, Functions::ZeroFunction<dim>(), constraints);
  constraints.close();

  deallog << "Testing Portable::MatrixFree on level " << test_level
          << std::endl;
  deallog << "Number of dofs on level " << test_level << ": "
          << dof_handler.n_dofs(test_level) << std::endl;

  // Initialize Portable::MatrixFree for level cells
  Portable::MatrixFree<dim, Number> mf_level_data;
  {
    typename Portable::MatrixFree<dim, Number>::AdditionalData additional_data;
    additional_data.mapping_update_flags = update_gradients | update_JxW_values;
    additional_data.mg_level             = test_level;
    const QGauss<1> quad(n_q_points_1d);

    mf_level_data.reinit(
      mapping, dof_handler, constraints, quad, additional_data);
  }

  // Initialize Portable::MatrixFree for active cells (for comparison)
  Portable::MatrixFree<dim, Number> mf_data;
  {
    typename Portable::MatrixFree<dim, Number>::AdditionalData additional_data;
    additional_data.mapping_update_flags = update_gradients | update_JxW_values;
    additional_data.mg_level             = numbers::invalid_unsigned_int;
    const QGauss<1> quad(n_q_points_1d);
    mf_data.reinit(mapping, dof_handler, constraints, quad, additional_data);
  }

  // Check that mg_level is stored correctly
  AssertThrow(mf_level_data.get_mg_level() == test_level,
              ExcMessage("mg_level not correctly stored"));

  // Create test vectors
  const unsigned int n_dofs = dof_handler.n_dofs(test_level);

  // Check that number of dofs match
  AssertThrow(n_dofs == dof_handler.n_dofs(),
              ExcMessage("Number of dofs do not match!"));

  LinearAlgebra::distributed::Vector<Number, MemorySpace::Default> src(n_dofs);
  LinearAlgebra::distributed::Vector<Number, MemorySpace::Default> dst_device(
    n_dofs);
  LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>
    dst_device_ref(n_dofs);

  // Fill source vector with random values
  LinearAlgebra::ReadWriteVector<Number> src_rw(n_dofs);
  for (unsigned int i = 0; i < n_dofs; ++i)
    {
      if (constraints.is_constrained(i))
        continue;
      const double entry = Testing::rand() / (double)RAND_MAX;
      src_rw(i)          = entry;
    }

  src.import_elements(src_rw, VectorOperation::insert);

  // Apply level operator using Portable::MatrixFree
  mf_level_data.cell_loop(
    LaplaceOperator<dim, fe_degree, Number, n_q_points_1d>(), src, dst_device);
  Kokkos::fence();

  // Apply active operator using Portable::MatrixFree for reference
  mf_data.cell_loop(LaplaceOperator<dim, fe_degree, Number, n_q_points_1d>(),
                    src,
                    dst_device_ref);
  Kokkos::fence();
  // Compare results

  LinearAlgebra::distributed::Vector<Number, MemorySpace::Default> diff(n_dofs);
  diff = dst_device;
  diff -= dst_device_ref;

  Number error_norm = 0.0;

  deallog << "Norm of difference: " << diff.l2_norm() << std::endl << std::endl;
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
