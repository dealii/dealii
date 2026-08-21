// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2019 by the deal.II authors
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

#include <deal.II/base/convergence_table.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/assemble_operator.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include <chrono>

using namespace dealii;

template <int dim,
          int fe_degree,
          int n_points,
          int n_components,
          typename Number>
class LinearOpertor
{
public:
  LinearOpertor(const MatrixFree<dim, Number> &matrix_free)
    : matrix_free(matrix_free)
  {}

  template <typename VectorType>
  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    FEEvaluation<dim, fe_degree, n_points, n_components, Number> phi(
      matrix_free);

    matrix_free.template cell_loop<VectorType, VectorType>(
      [&](const auto &, auto &dst, const auto &src, const auto cells) {
        for (auto cell = cells.first; cell < cells.second; cell++)
          {
            phi.reinit(cell);
            phi.gather_evaluate(src, false, true);
            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              phi.submit_gradient(phi.get_gradient(q), q);
            phi.integrate_scatter(false, true, dst);
          }
      },
      dst,
      src);
  }

  const MatrixFree<dim, Number> &matrix_free;
};

template <int dim,
          int fe_degree,
          int n_points,
          int n_components,
          typename Number>
void
do_test(ConvergenceTable &table, const unsigned int n_refinements)
{
  using VectorType = Vector<Number>;

  // 1) setup triangulation
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(n_refinements);

  // 2) create dof-handler
  const FE_Q<dim>     fe_q(fe_degree);
  const FESystem<dim> fe(fe_q, n_components);

  // setup dof-handlers
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  const unsigned int n_dofs = dof_handler.n_dofs();

  // 3) create constraints
  AffineConstraints<Number> constraint;
  if (false)
    {
      // hanging nodes
      DoFTools::make_hanging_node_constraints(dof_handler, constraint);

      // Dirichlet boundary
      VectorTools::interpolate_boundary_values(
        dof_handler, 0, Functions::ZeroFunction<dim>(n_components), constraint);
    }
  constraint.close();

  // 4) create matrix free
  typename MatrixFree<dim, Number>::AdditionalData additional_data;
  additional_data.tasks_parallel_scheme =
    decltype(additional_data)::TasksParallelScheme::none;
  additional_data.mapping_update_flags = update_values | update_gradients;

  MappingQ<dim> mapping(1);
  QGauss<1>     quad(n_points);

  MatrixFree<dim, Number> matrix_free;
  matrix_free.reinit(mapping, dof_handler, constraint, quad, additional_data);

  // 5) create linear operator
  LinearOpertor<dim, fe_degree, n_points, n_components, Number> linear_operator(
    matrix_free);

  // 6) create utility objects
  std::shared_ptr<MatrixFreeTools::GraphCache> cache;

  // ... vectors
  VectorType register_1(n_dofs);
  VectorType register_2(n_dofs);

  // .. pattern
  DynamicSparsityPattern dynamic_sparsity_pattern(n_dofs, n_dofs);
  DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern);
  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(dynamic_sparsity_pattern);

  // 7) tests

  // full matrix
  {
    auto               start = std::chrono::steady_clock::now();
    FullMatrix<Number> matrix(n_dofs, n_dofs);
    MatrixFreeTools::assemble_operator(
      matrix, linear_operator, sparsity_pattern, register_1, register_2, cache);
    auto                          end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    table.add_value("t_full", elapsed_seconds.count());
    table.set_scientific("t_full", true);
  }

  // sparse matrix
  {
    auto                 start = std::chrono::steady_clock::now();
    SparseMatrix<Number> matrix(sparsity_pattern);
    MatrixFreeTools::assemble_operator(
      matrix, linear_operator, sparsity_pattern, register_1, register_2, cache);
    auto                          end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    table.add_value("t_sm", elapsed_seconds.count());
    table.set_scientific("t_sm", true);
  }

  // reference
  {
    auto start = std::chrono::steady_clock::now();
    for (unsigned int i = 0; i < cache->num_colors; i++)
      linear_operator.vmult(register_2, register_1);
    auto                          end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    table.add_value("t_vmult", elapsed_seconds.count());
    table.set_scientific("t_vmult", true);
  }

  table.add_value("dim", dim);
  table.add_value("degree", fe_degree);
  table.add_value("n_dofs", n_dofs);
  table.add_value("n_vmults", cache->num_colors);
  table.add_value("ratio", 1.0 * n_dofs / cache->num_colors);

  // matrix.print_formatted (std::cout);
}

template <int dim, typename Number>
void
do_test_dim(ConvergenceTable &table)
{
  for (unsigned int i = 0; i < 4; i++)
    {
      do_test<dim, 1, 2, 1, Number>(table, i);   // linear (scalar)
      do_test<dim, 2, 3, 1, Number>(table, i);   // quadratic (scalar)
      do_test<dim, 1, 2, dim, Number>(table, i); // linear (vector)
      do_test<dim, 2, 3, dim, Number>(table, i); // quadratic (vector)
    }
}

int
main()
{
  ConvergenceTable table;

  do_test_dim<2, double>(table); // 2D
  do_test_dim<3, double>(table); // 3D

  table.write_text(std::cout);
}
