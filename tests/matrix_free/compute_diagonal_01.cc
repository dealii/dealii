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


// Test MatrixFreeTools::compute_diagonal() for a Laplace operator
// If the flag 'use_categorization = true' is set, additional categorization
// of cells according to their material_id is performed and enforced in
// MatrixFree via MF::AD::cell_vectorization_categories_strict See
// https://github.com/dealii/dealii/issues/16250 for more details.

#include "compute_diagonal_util.h"

template <int dim,
          int fe_degree,
          int n_points                 = fe_degree + 1,
          int n_components             = 1,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
void
test(bool use_categorization = false)
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tria, -numbers::PI / 2, numbers::PI / 2);

  tria.refine_global();
  for (auto &cell : tria.active_cell_iterators())
    if (cell->is_active() && cell->is_locally_owned() &&
        cell->center()[0] < 0.0)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  if (use_categorization)
    {
      for (auto &cell : tria.active_cell_iterators())
        if (cell->is_active() && cell->is_locally_owned() &&
            cell->center()[0] < 0.0)
          cell->set_material_id(42);
        else
          cell->set_material_id(0);
    }
  const FE_Q<dim>     fe_q(fe_degree);
  const FESystem<dim> fe(fe_q, n_components);

  // setup dof-handlers
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<Number> constraint;
  DoFTools::make_hanging_node_constraints(dof_handler, constraint);

  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           Functions::ZeroFunction<dim, Number>(
                                             n_components),
                                           constraint);

  constraint.close();

  typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
    additional_data;
  additional_data.mapping_update_flags = update_values | update_gradients;

  if (use_categorization)
    {
      additional_data.cell_vectorization_category.resize(tria.n_active_cells());
      for (const auto &cell : tria.active_cell_iterators())
        if (cell->is_locally_owned())
          additional_data
            .cell_vectorization_category[cell->active_cell_index()] =
            cell->material_id();

      additional_data.cell_vectorization_categories_strict = true;
    }

  MappingQ<dim> mapping(1);
  QGauss<1>     quad(fe_degree + 1);

  MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
  matrix_free.reinit(mapping, dof_handler, constraint, quad, additional_data);


  Test<dim, fe_degree, n_points, n_components, Number, VectorizedArrayType>
    test(matrix_free,
         constraint,
         [](FEEvaluation<dim,
                         fe_degree,
                         n_points,
                         n_components,
                         Number,
                         VectorizedArrayType> &phi) {
           phi.evaluate(EvaluationFlags::gradients);
           for (unsigned int q = 0; q < phi.n_q_points; ++q)
             phi.submit_gradient(phi.get_gradient(q), q);
           phi.integrate(EvaluationFlags::gradients);
         });

  test.do_test();
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test<2, 1, 2, 1>(); // scalar
  test<2, 1, 2, 2>(); // vector

  test<2, 1, 2, 1, float>(); // scalar
  test<2, 1, 2, 2, float>(); // vector

  // Same as above, but testing additional categorization of vectorized cell
  // batches
  test<2, 1, 2, 1, float>(true); // scalar
  test<2, 1, 2, 2, float>(true); // vector
}
