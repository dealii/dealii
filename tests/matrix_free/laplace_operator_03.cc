// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// the same as laplace_operator_01 (excluding the extra detection tests), but
// heterogeneous Laplace operator with a single constant coefficient per cell

#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <deal.II/matrix_free/operators.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"


template <int dim, int fe_degree>
void
test()
{
  using number = double;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);
  for (const auto &cell :
       tria.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
    if (cell->center().norm() < 0.2)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  if (dim < 3 && fe_degree < 2)
    tria.refine_global(2);
  else
    tria.refine_global(1);
  if (tria.begin(tria.n_levels() - 1)->is_locally_owned())
    tria.begin(tria.n_levels() - 1)->set_refine_flag();
  if (tria.last()->is_locally_owned())
    tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  for (unsigned int i = 0; i < 10 - 3 * dim; ++i)
    {
      unsigned int counter = 0;
      for (const auto &cell :
           tria.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
        if (counter++ % (7 - i) == 0)
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);

  const IndexSet &owned_set    = dof.locally_owned_dofs();
  const IndexSet  relevant_set = DoFTools::extract_locally_relevant_dofs(dof);

  AffineConstraints<double> constraints(owned_set, relevant_set);
  DoFTools::make_hanging_node_constraints(dof, constraints);
  VectorTools::interpolate_boundary_values(dof,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  constraints.close();

  deallog << "Testing " << dof.get_fe().get_name() << std::endl;
  // std::cout << "Number of cells: " << tria.n_global_active_cells() <<
  // std::endl; std::cout << "Number of degrees of freedom: " << dof.n_dofs() <<
  // std::endl; std::cout << "Number of constraints: " <<
  // constraints.n_constraints() << std::endl;

  std::shared_ptr<MatrixFree<dim, number>> mf_data(
    new MatrixFree<dim, number>());
  {
    const QGauss<1>                                  quad(fe_degree + 1);
    typename MatrixFree<dim, number>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim, number>::AdditionalData::none;
    data.tasks_block_size      = 7;
    data.mapping_update_flags =
      update_quadrature_points | update_gradients | update_JxW_values;
    mf_data->reinit(MappingQ1<dim>{}, dof, constraints, quad, data);
  }

  std::shared_ptr<Table<2, VectorizedArray<number>>> coefficient;
  coefficient = std::make_shared<Table<2, VectorizedArray<number>>>();
  {
    coefficient->reinit(mf_data->n_cell_batches(), 1);
    for (unsigned int cell = 0; cell < mf_data->n_cell_batches(); ++cell)
      for (unsigned int v = 0;
           v < mf_data->n_active_entries_per_cell_batch(cell);
           ++v)
        (*coefficient)(cell, 0)[v] =
          1. + std::sqrt(static_cast<double>(
                 mf_data->get_cell_iterator(cell, v)->active_cell_index()));
  }

  MatrixFreeOperators::LaplaceOperator<
    dim,
    fe_degree,
    fe_degree + 1,
    1,
    LinearAlgebra::distributed::Vector<number>>
    mf;
  mf.initialize(mf_data);
  mf.set_coefficient(coefficient);
  mf.compute_diagonal();
  LinearAlgebra::distributed::Vector<number> in, out, ref;
  mf_data->initialize_dof_vector(in);
  out.reinit(in);
  ref.reinit(in);

  for (unsigned int i = 0; i < in.locally_owned_size(); ++i)
    {
      const unsigned int glob_index = owned_set.nth_index_in_set(i);
      if (constraints.is_constrained(glob_index))
        continue;
      in.local_element(i) = random_value<double>();
    }

  mf.vmult(out, in);


  // assemble trilinos sparse matrix with
  // (v, u) for reference
  TrilinosWrappers::SparseMatrix sparse_matrix;
  {
    TrilinosWrappers::SparsityPattern csp(owned_set, MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof,
                                    csp,
                                    constraints,
                                    true,
                                    Utilities::MPI::this_mpi_process(
                                      MPI_COMM_WORLD));
    csp.compress();
    sparse_matrix.reinit(csp);
  }
  {
    QGauss<dim> quadrature_formula(fe_degree + 1);

    FEValues<dim> fe_values(dof.get_fe(),
                            quadrature_formula,
                            update_gradients | update_JxW_values |
                              update_quadrature_points);

    const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(),
                                                   endc = dof.end();
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
        {
          cell_matrix = 0;
          fe_values.reinit(cell);

          const double coefficient =
            1. + std::sqrt(static_cast<double>(cell->active_cell_index()));

          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  cell_matrix(i, j) += (fe_values.shape_grad(i, q_point) *
                                        fe_values.shape_grad(j, q_point)) *
                                       coefficient * fe_values.JxW(q_point);
              }

          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_matrix,
                                                 local_dof_indices,
                                                 sparse_matrix);
        }
  }
  sparse_matrix.compress(VectorOperation::add);

  sparse_matrix.vmult(ref, in);
  out -= ref;
  const double diff_norm = out.linfty_norm();

  deallog << "Norm of difference: " << diff_norm << std::endl << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  mpi_initlog();

  deallog << std::setprecision(4);

  deallog.push("2d");
  test<2, 1>();
  test<2, 2>();
  deallog.pop();

  deallog.push("3d");
  test<3, 1>();
  test<3, 2>();
  deallog.pop();
}
