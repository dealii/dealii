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


// Same as compute_diagonal_01 but with two different DoFHandlers attached to
// MatrixFree, one of which is in hp-mode

#include "compute_diagonal_util.h"

template <int dim,
          int fe_degree                = 1,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tria, -numbers::PI / 2, numbers::PI / 2);

  tria.refine_global();
  for (auto &cell : tria.active_cell_iterators())
    if (cell->is_active() && cell->is_locally_owned() &&
        cell->center()[0] < 0.0)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FE_Q<dim> fe0(fe_degree);
  FE_Q<dim> fe1(fe_degree + 1);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(fe0);
  fe_collection.push_back(fe1);

  DoFHandler<dim> dof0(tria);
  dof0.distribute_dofs(fe0);

  DoFHandler<dim> dof1(tria);
  unsigned int    counter = 0;
  for (const auto &cell : dof1.active_cell_iterators())
    cell->set_active_fe_index(++counter % 2);
  dof1.distribute_dofs(fe_collection);

  std::vector<const DoFHandler<dim> *> dof(2);
  dof[0] = &dof0;
  dof[1] = &dof1;

  AffineConstraints<Number> constraint0;
  DoFTools::make_hanging_node_constraints(*dof[0], constraint0);
  VectorTools::interpolate_boundary_values(
    *dof[0], 0, Functions::ZeroFunction<dim, Number>(), constraint0);
  constraint0.close();

  AffineConstraints<Number> constraint1;
  DoFTools::make_hanging_node_constraints(*dof[1], constraint1);
  VectorTools::interpolate_boundary_values(
    *dof[1], 0, Functions::ZeroFunction<dim, Number>(), constraint1);
  constraint1.close();

  std::vector<const AffineConstraints<Number> *> constraints(2);
  constraints[0] = &constraint0;
  constraints[1] = &constraint1;


  typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
    additional_data;
  additional_data.mapping_update_flags = update_values | update_gradients;

  MappingQ<dim> mapping(1);

  std::vector<Quadrature<dim>> quad;
  for (unsigned int no = 0; no < 2; ++no)
    quad.push_back(QGauss<dim>(fe_degree + 1 + no));

  MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
  matrix_free.reinit(mapping, dof, constraints, quad, additional_data);


  Test<dim, -1, 0, 1, Number, VectorizedArrayType> test(
    matrix_free,
    *constraints[0],
    [](FEEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType> &phi) {
      phi.evaluate(EvaluationFlags::gradients);
      for (unsigned int q = 0; q < phi.n_q_points; ++q)
        phi.submit_gradient(phi.get_gradient(q), q);
      phi.integrate(EvaluationFlags::gradients);
    },
    0,
    0);

  test.do_test();
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test<2>(); // scalar
}
