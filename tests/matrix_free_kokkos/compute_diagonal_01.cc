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


// Test MatrixFreeTools::compute_diagonal() for a Laplace operator

#include "compute_diagonal_util.h"

template <int dim,
          int fe_degree,
          int n_points     = fe_degree + 1,
          int n_components = 1,
          typename Number  = double>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tria, -numbers::PI / 2, numbers::PI / 2);

  tria.refine_global(2);
  for (auto &cell : tria.active_cell_iterators())
    if (cell->is_active() && cell->is_locally_owned() &&
        cell->center()[0] < 0.0)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

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

  typename Portable::MatrixFree<dim, Number>::AdditionalData additional_data;
  additional_data.mapping_update_flags = update_values | update_gradients;

  MappingQ<dim> mapping(1);
  QGauss<1>     quad(fe_degree + 1);

  Portable::MatrixFree<dim, Number> matrix_free;
  matrix_free.reinit(mapping, dof_handler, constraint, quad, additional_data);


  Test<dim, fe_degree, n_points, n_components, Number> test(matrix_free,
                                                            constraint);

  deallog << "dim: " << dim << ", fe_degree: " << fe_degree
          << ", n_points: " << n_points << ", n_components: " << n_components
          << ", Number: "
          << (std::is_same_v<Number, double> ? "double" : "float") << std::endl;

  test.do_test();
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test<2, 1, 2, 1>();        // scalar
  test<2, 1, 2, 1, float>(); // scalar
  test<3, 1, 2, 1>();        // scalar
  test<3, 1, 2, 1, float>(); // scalar

  test<2, 1, 2, 2>();        // vectorial
  test<2, 1, 2, 2, float>(); // vectorial
  test<3, 1, 2, 2>();        // vectorial
  test<3, 1, 2, 2, float>(); // vectorial
}
