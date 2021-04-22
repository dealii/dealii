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


// Same as compute_diagonal_01 but for hp.

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
    if (cell->active() && cell->is_locally_owned() && cell->center()[0] < 0.0)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();


  hp::FECollection<dim> fes{FE_Q<dim>(1), FE_Q<dim>(2)};

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


  Test<dim, -1, 0, n_components, Number, VectorizedArrayType> test(
    matrix_free,
    constraint,
    [](FEEvaluation<dim, -1, 0, n_components, Number, VectorizedArrayType>
         &phi) {
      phi.evaluate(false, true);
      for (unsigned int q = 0; q < phi.n_q_points; ++q)
        phi.submit_gradient(phi.get_gradient(q), q);
      phi.integrate(false, true);
    });

  test.do_test();
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test<2>();
}
