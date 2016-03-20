// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// DoFTools::locally_owned_dofs_per_subdomain wrongly sized the array
// it returns when we have some processors that don't own any cells

#include "../tests.h"
#include <deal.II/grid/grid_generator.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>


template <int dim>
void test ()
{
  parallel::shared::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(triangulation, -1.0, 1.0);

  FE_Q<dim> fe(1);
  DoFHandler<dim> dof_handler (triangulation);
  dof_handler.distribute_dofs(fe);

  // this used to crash here:
  const IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();

  deallog << "dim=" << dim
          << std::endl
          << "rank=" << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
          << std::endl
          << "  n_dofs=" << dof_handler.n_dofs()
          << std::endl
          << "  n_locally_owned_dofs="
          << locally_owned_dofs.n_elements()
          << std::endl;
}

int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, testing_max_num_threads());

  MPILogInitAll all;

  test<2>();
  test<3>();
}
