// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2024 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Test MatrixFreeTools::compute_matrix() for Portable::MatrixFree against
// FEValues assembly of a Laplace-plus-mass operator. Cover scalar and coupled
// components, homogeneous Dirichlet and hanging-node constraints, deformed
// meshes, coloring, and repeated assembly into the same matrix.

#include <deal.II/base/function_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>
#include <deal.II/lac/trilinos_tpetra_sparsity_pattern.h>

#include <deal.II/matrix_free/portable_fe_evaluation.h>
#include <deal.II/matrix_free/portable_matrix_free.h>
#include <deal.II/matrix_free/tools.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"


template <int dim, int degree, int n_components>
struct CellOperation
{
  DEAL_II_HOST_DEVICE void
  operator()(Portable::FEEvaluation<dim, degree, degree + 1, n_components> *phi,
             const int q) const
  {
    auto value = phi->get_value(q);
    if constexpr (n_components > 1)
      {
        const auto original = value;
        for (unsigned int c = 0; c < n_components; ++c)
          value[c] += 0.25 * original[(c + 1) % n_components];
      }
    phi->submit_value(value, q);
    phi->submit_gradient(phi->get_gradient(q), q);
  }
};


template <int dim, int degree, int n_components = 1>
void
test(const bool use_coloring, const bool deform, const bool adaptive = false)
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1);
  if (adaptive)
    {
      triangulation.begin_active()->set_refine_flag();
      triangulation.execute_coarsening_and_refinement();
    }
  if (deform)
    GridTools::transform(
      [](const Point<dim> &p) {
        Point<dim> result = p;
        result[0] += 0.1 * p[0] * p[1];
        return result;
      },
      triangulation);

  FESystem<dim>   fe(FE_Q<dim>{degree}, n_components);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  auto constraints = DoFTools::make_default_constraints(dof_handler);
  VectorTools::interpolate_boundary_values(
    dof_handler, 0, Functions::ZeroFunction<dim>(n_components), constraints);
  constraints.close();

  DynamicSparsityPattern dynamic_sparsity(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dynamic_sparsity,
                                  constraints,
                                  false);
  SparsityPattern sparsity;
  sparsity.copy_from(dynamic_sparsity);

  LinearAlgebra::TpetraWrappers::SparseMatrix<double, MemorySpace::Default>
    matrix;
  matrix.reinit(sparsity);

  Portable::MatrixFree<dim, double>                          matrix_free;
  typename Portable::MatrixFree<dim, double>::AdditionalData additional_data;
  additional_data.mapping_update_flags = update_gradients | update_JxW_values;
  additional_data.use_coloring         = use_coloring;
  matrix_free.reinit(dof_handler,
                     constraints,
                     QGauss<1>(degree + 1),
                     additional_data);

  // Repeated assembly must clear the previous entries, including constrained
  // diagonals, and release all device views before returning.
  for (unsigned int repetition = 0; repetition < 2; ++repetition)
    MatrixFreeTools::compute_matrix<dim, degree, degree + 1, n_components>(
      matrix_free,
      matrix,
      CellOperation<dim, degree, n_components>(),
      EvaluationFlags::values | EvaluationFlags::gradients,
      EvaluationFlags::values | EvaluationFlags::gradients);

  // Compare every assembled entry against independent host FEValues assembly.
  SparseMatrix<double> reference(sparsity);
  FEValues<dim> fe_values(fe,
                          QGauss<dim>(degree + 1),
                          update_values | update_gradients | update_JxW_values);
  FullMatrix<double> cell_matrix(fe.n_dofs_per_cell(), fe.n_dofs_per_cell());
  std::vector<types::global_dof_index> indices(fe.n_dofs_per_cell());
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      cell_matrix = 0.;
      for (const unsigned int q : fe_values.quadrature_point_indices())
        for (const unsigned int i : fe_values.dof_indices())
          for (const unsigned int j : fe_values.dof_indices())
            for (unsigned int c = 0; c < n_components; ++c)
              {
                cell_matrix(i, j) +=
                  (fe_values.shape_value_component(i, q, c) *
                     fe_values.shape_value_component(j, q, c) +
                   fe_values.shape_grad_component(i, q, c) *
                     fe_values.shape_grad_component(j, q, c)) *
                  fe_values.JxW(q);
                if constexpr (n_components > 1)
                  cell_matrix(i, j) +=
                    0.25 * fe_values.shape_value_component(i, q, c) *
                    fe_values.shape_value_component(j,
                                                    q,
                                                    (c + 1) % n_components) *
                    fe_values.JxW(q);
              }
      cell->get_dof_indices(indices);
      constraints.distribute_local_to_global(cell_matrix, indices, reference);
    }
  for (const auto &line : constraints.get_lines())
    reference.set(line.index, line.index, 1.);

  const auto host_matrix = matrix.trilinos_matrix().getLocalMatrixHost();
  double     error       = 0.;
  for (unsigned int row = 0; row < dof_handler.n_dofs(); ++row)
    for (auto entry = host_matrix.graph.row_map(row);
         entry < host_matrix.graph.row_map(row + 1);
         ++entry)
      error =
        std::max(error,
                 std::abs(host_matrix.values(entry) -
                          reference.el(row, host_matrix.graph.entries(entry))));
  AssertThrow(error < 1.e-11,
              ExcMessage("Device assembly differs from FEValues"));
  deallog << "dim=" << dim << " degree=" << degree
          << " coloring=" << use_coloring << " deform=" << deform;
  if (n_components > 1)
    deallog << " components=" << n_components;
  if (adaptive)
    deallog << " adaptive=" << adaptive;
  deallog << ": OK" << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());
  initlog();

  for (const bool use_coloring : {false, true})
    {
      test<2, 1>(use_coloring, false);
      test<2, 2>(use_coloring, true);
      test<3, 1>(use_coloring, false);
      test<3, 2>(use_coloring, true);
      test<2, 2, 2>(use_coloring, true, true);
      test<2, 2>(use_coloring, true, true);
    }
}
