// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Test MatrixFreeTools::compute_diagonal() for a System with 2 DoFHandlers

#include "compute_diagonal_util.h"

template <int dim, int fe_degree, int n_points, typename Number = double>
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

  const unsigned int degree_u = fe_degree + 1;
  const unsigned int degree_p = fe_degree;

  FESystem<dim>   fe_u(FE_Q<dim>(degree_u), dim);
  FE_Q<dim>       fe_p(degree_p);
  DoFHandler<dim> dof_u(tria);
  DoFHandler<dim> dof_p(tria);
  dof_u.distribute_dofs(fe_u);
  dof_p.distribute_dofs(fe_p);

  const IndexSet &owned_set_u = dof_u.locally_owned_dofs();
  const IndexSet  relevant_set_u =
    DoFTools::extract_locally_relevant_dofs(dof_u);
  AffineConstraints<Number> constraints_u(owned_set_u, relevant_set_u);
  DoFTools::make_hanging_node_constraints(dof_u, constraints_u);
  VectorTools::interpolate_boundary_values(
    dof_u, 0, Functions::ZeroFunction<dim, Number>(dim), constraints_u);
  constraints_u.close();

  const IndexSet &owned_set_p = dof_p.locally_owned_dofs();
  const IndexSet  relevant_set_p =
    DoFTools::extract_locally_relevant_dofs(dof_p);
  AffineConstraints<Number> constraints_p(owned_set_p, relevant_set_p);
  DoFTools::make_hanging_node_constraints(dof_p, constraints_p);
  constraints_p.close();

  std::vector<const DoFHandler<dim> *> dof_handlers          = {&dof_u, &dof_p};
  std::vector<const AffineConstraints<Number> *> constraints = {&constraints_u,
                                                                &constraints_p};


  typename Portable::MatrixFree<dim, Number>::AdditionalData additional_data;
  additional_data.mapping_update_flags = update_values | update_gradients;

  MappingQ<dim> mapping(1);
  QGauss<1>     quad(fe_degree + 2);

  Portable::MatrixFree<dim, Number> matrix_free;
  matrix_free.reinit(mapping, dof_handlers, constraints, quad, additional_data);


  {
    Test<dim, degree_u, n_points, dim, Number> test(matrix_free,
                                                    constraints_u,
                                                    0 /* velocity */);

    deallog << "dim: " << dim << ", fe_degree: " << fe_degree
            << ", n_points: " << n_points << ", n_components: " << dim
            << ", dof_handler_index: " << 0 << ", Number: "
            << (std::is_same_v<Number, double> ? "double" : "float")
            << std::endl;

    test.do_test();
  }
  {
    Test<dim, degree_p, n_points, 1, Number> test(matrix_free,
                                                  constraints_p,
                                                  1 /* pressure */);

    deallog << "dim: " << dim << ", fe_degree: " << fe_degree
            << ", n_points: " << n_points << ", n_components: " << 1
            << ", dof_handler_index: " << 1 << ", Number: "
            << (std::is_same_v<Number, double> ? "double" : "float")
            << std::endl;

    test.do_test();
  }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  // Q2-Q1:
  test<2, 1, 3>();
  test<2, 1, 3, float>();
  test<3, 1, 3>();
  test<3, 1, 3, float>();
}
