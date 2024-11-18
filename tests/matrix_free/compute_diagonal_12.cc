// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Based on compute_diagonal_04 but uses FE_Nothing

#include "compute_diagonal_util.h"

template <int dim,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
void
test()
{
  const unsigned int n_components = 1;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tria, -numbers::PI / 2, numbers::PI / 2);

  tria.refine_global();
  for (auto &cell : tria.active_cell_iterators())
    if (cell->is_active() && cell->is_locally_owned() &&
        cell->center()[0] < 0.0)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();


  hp::FECollection<dim> fes{FE_Q<dim>(1), FE_Nothing<dim>(1)};

  // setup dof-handlers
  DoFHandler<dim> dof_handler(tria);

  unsigned int counter = 0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    cell->set_active_fe_index(++counter % 2);

  dof_handler.distribute_dofs(fes);

  AffineConstraints<Number> constraint;
  DoFTools::make_hanging_node_constraints(dof_handler, constraint);

  VectorTools::interpolate_boundary_values(
    dof_handler, 0, Functions::ZeroFunction<dim>(n_components), constraint);

  constraint.close();

  typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
    additional_data;
  additional_data.mapping_update_flags = update_values | update_gradients;

  MappingQ<dim> mapping(1);
  QGauss<1>     quad(2);

  MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
  matrix_free.reinit(mapping, dof_handler, constraint, quad, additional_data);

  LinearAlgebra::distributed::Vector<Number> diagonal;
  matrix_free.initialize_dof_vector(diagonal);

  dealii::MatrixFreeTools::compute_diagonal<
    dim,
    -1,
    0,
    n_components,
    Number,
    VectorizedArrayType,
    LinearAlgebra::distributed::Vector<Number>>(
    matrix_free,
    diagonal,
    [](FEEvaluation<dim, -1, 0, n_components, Number, VectorizedArrayType>
         &phi) {
      phi.evaluate(EvaluationFlags::gradients);
      for (unsigned int q = 0; q < phi.n_q_points; ++q)
        phi.submit_gradient(phi.get_gradient(q), q);
      phi.integrate(EvaluationFlags::gradients);
    });

  diagonal.print(deallog.get_file_stream());
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test<2>();
}
