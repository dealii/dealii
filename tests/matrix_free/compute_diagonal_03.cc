// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2021 by the deal.II authors
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



// Test MatrixFreeTools::compute_diagonal() for scalar and vector Laplace
// operator with variable coefficients.

#include "compute_diagonal_util.h"

template <int dim,
          int fe_degree,
          int n_points                 = fe_degree + 1,
          int n_components             = 1,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);

  GridGenerator::hyper_ball(tria);
  tria.refine_global(1);

  const FE_Q<dim>     fe_q(fe_degree);
  const FESystem<dim> fe(fe_q, n_components);

  // setup dof-handlers
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<Number> constraint;
  DoFTools::make_hanging_node_constraints(dof_handler, constraint);

  VectorTools::interpolate_boundary_values(
    dof_handler, 0, Functions::ZeroFunction<dim>(n_components), constraint);

  constraint.close();

  typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
    additional_data;
  additional_data.mapping_update_flags =
    update_values | update_gradients | update_quadrature_points;

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
           phi.evaluate(false, true);
           for (unsigned int q = 0; q < phi.n_q_points; ++q)
             {
               auto                quadrature_point = phi.quadrature_point(q);
               VectorizedArrayType coefficient;

               for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
                 {
                   Point<dim, Number> point;
                   for (unsigned int d = 0; d < dim; ++d)
                     point[d] = quadrature_point[d][v];
                   coefficient[v] = (point.square() < 0.5 * 0.5) ? 20 : 1;
                 }

               phi.submit_gradient(coefficient * phi.get_gradient(q), q);
             }
           phi.integrate(false, true);
         });

  test.do_test();
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test<2, 1, 2, 1>(); // scalar
  test<2, 2, 3, 1>(); // scalar
  test<2, 3, 4, 1>(); // scalar
  test<3, 1, 2, 1>(); // vector
  test<3, 2, 3, 1>(); // vector
  test<3, 3, 4, 1>(); // vector
}
