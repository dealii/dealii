// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2019 by the deal.II authors
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



// test that DoFHandler::n_boundary_dofs() yields the correct
// results for parallel::distributed::Triangulations

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(3);

  FE_Q<dim>       fe(2);
  DoFHandler<dim> dof_handler(triangulation);

  dof_handler.distribute_dofs(fe);

  const unsigned int local_boundary_dofs = dof_handler.n_boundary_dofs();
  const unsigned int global_boundary_dofs =
    Utilities::MPI::sum(local_boundary_dofs, MPI_COMM_WORLD);
  deallog << global_boundary_dofs << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
