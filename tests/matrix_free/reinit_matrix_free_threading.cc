// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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

// Verify that MatrixFree::reinit() works without errors if run with
// multithreadingf and multiple ranks.


#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/matrix_free/matrix_free.h>

#include <iostream>

#include "../tests.h"


template <int dim, typename Number>
void
do_test(const unsigned int degree, const unsigned int n_refinements)
{
  deallog << "Create Triangulation with " << n_refinements << " refinements..."
          << std::endl;
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(n_refinements);


  deallog << "Create DoFHandler..." << std::endl;
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FE_DGQ<dim>(degree));

  deallog << "Setup AffineConstraints..." << std::endl;
  AffineConstraints<Number> constraints;
  constraints.close();

  deallog << "Set MatrixFree AdditionalData..." << std::endl;
  typename MatrixFree<dim, Number>::AdditionalData data;
  data.mapping_update_flags_inner_faces = update_values;

  deallog << "Reinit MatrixFree..." << std::endl;
  MatrixFree<dim, Number> matrix_free;
  matrix_free.reinit(
    MappingQ1<dim>(), dof_handler, constraints, QGauss<dim>(degree + 1), data);

  deallog << "OK" << std::endl << std::endl;
}



int
main(int argc, char **argv)
{
  // Init with threading enabled
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  mpi_initlog();

  do_test<2, double>(3, 3);
  do_test<3, double>(2, 2);

  deallog << std::endl;
}
