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


// Test MatrixFreeTools::compute_diagonal() for a Laplace operator on a simplex
// mesh.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>

#include "compute_diagonal_util.h"

template <int dim,
          int n_components             = 1,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
void
test()
{
  Triangulation<dim> tria;

  GridGenerator::subdivided_hyper_cube_with_simplices(tria, 2);

  const unsigned int fe_degree = 2;
  const unsigned int n_points  = 3;

  const FE_SimplexP<dim> fe_q(fe_degree);
  const FESystem<dim>    fe(fe_q, n_components);

  // setup dof-handlers
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<Number> constraint;

  MappingFE<dim> mapping(FE_SimplexP<dim>{1});

  VectorTools::interpolate_boundary_values(mapping,
                                           dof_handler,
                                           0,
                                           Functions::ZeroFunction<dim>(
                                             n_components),
                                           constraint);

  constraint.close();

  typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
    additional_data;
  additional_data.mapping_update_flags = update_values | update_gradients;
  QGaussSimplex<dim> quad(fe_degree + 1);

  MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
  matrix_free.reinit(mapping, dof_handler, constraint, quad, additional_data);


  Test<dim, -1, 0, n_components, Number, VectorizedArrayType> test(
    matrix_free,
    constraint,
    [](FEEvaluation<dim, -1, 0, n_components, Number, VectorizedArrayType>
         &phi) {
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

  test<2, 1>(); // scalar
  test<2, 2>(); // vector
}
