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



// Test MatrixFreeTools::compute_diagonal() for a vector Laplace operator
// and compute_no_normal_flux_constraints.

#include "compute_diagonal_util.h"

template <int dim>
class BoundaryValues : public Function<dim>
{
public:
  BoundaryValues()
    : Function<dim>(dim)
  {}
  virtual void
  vector_value(const Point<dim> &p, Vector<double> &value) const;
};
template <int dim>
void
BoundaryValues<dim>::vector_value(const Point<dim> &p,
                                  Vector<double> &  values) const
{
  (void)p;
  for (unsigned int i = 0; i < values.size(); ++i)
    values(i) = 0.0;
  return;
}

using namespace dealii;

template <int dim,
          int fe_degree,
          int n_points                 = fe_degree + 1,
          int n_components             = dim,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
void
test()
{
  Triangulation<dim> tria;

  GridGenerator::hyper_shell(tria, Point<dim>(), 1, 2, 0, true);

  tria.refine_global(0);


  const FE_Q<dim>     fe_q(fe_degree);
  const FESystem<dim> fe(fe_q, n_components);

  // setup dof-handlers
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<Number> constraints;
  IndexSet                  locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  constraints.reinit(locally_relevant_dofs);

  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

  std::set<unsigned int> dirichlet_boundary  = {0};
  std::set<unsigned int> tangential_boundary = {1};

  BoundaryValues<dim> boundary;
  for (auto bid : dirichlet_boundary)
    VectorTools::interpolate_boundary_values(dof_handler,
                                             bid,
                                             boundary,
                                             constraints);

  VectorTools::compute_no_normal_flux_constraints(dof_handler,
                                                  /* first_vector_component= */
                                                  0,
                                                  tangential_boundary,
                                                  constraints);
  constraints.close();

  typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
    additional_data;
  additional_data.mapping_update_flags = update_values | update_gradients;

  MappingQ<dim> mapping(1);
  QGauss<1>     quad(fe_degree + 1);

  MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
  matrix_free.reinit(mapping, dof_handler, constraints, quad, additional_data);


  Test<dim, fe_degree, n_points, n_components, Number, VectorizedArrayType>
    test(matrix_free,
         constraints,
         [](FEEvaluation<dim,
                         fe_degree,
                         n_points,
                         n_components,
                         Number,
                         VectorizedArrayType> &phi) {
           phi.evaluate(false, true, false);
           for (unsigned int q = 0; q < phi.n_q_points; ++q)
             {
               phi.submit_symmetric_gradient(2.0 *
                                               phi.get_symmetric_gradient(q),
                                             q);
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

  test<2, 1>();
}
