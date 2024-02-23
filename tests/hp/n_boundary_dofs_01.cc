// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test that DoFHandler::n_boundary_dofs() yields the correct
// results for parallel::distributed::Triangulations

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(3);

  FE_Q<dim>             fe(2);
  hp::FECollection<dim> fe_coll(fe);
  DoFHandler<dim>       dof_handler(triangulation);

  dof_handler.distribute_dofs(fe_coll);

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
