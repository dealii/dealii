// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that DoFHandler::clear() will be called
// whenever Triangulation::Signals::clear is triggered


#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


void
test_serial()
{
  constexpr const unsigned dim = 2;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);

  DoFHandler<dim> dh(tria);
  dh.distribute_dofs(FE_Q<dim>(1));

  tria.clear();
  deallog << "ndofs:" << dh.n_dofs() << std::endl;
}


void
test_parallel_shared()
{
  constexpr const unsigned dim = 2;

  parallel::shared::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);

  DoFHandler<dim> dh(tria);
  dh.distribute_dofs(FE_Q<dim>(1));

  tria.clear();
  deallog << "ndofs:" << dh.n_dofs() << std::endl;
}


void
test_parallel_distributed()
{
  constexpr const unsigned dim = 2;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);

  DoFHandler<dim> dh(tria);
  dh.distribute_dofs(FE_Q<dim>(1));

  tria.clear();
  deallog << "ndofs:" << dh.n_dofs() << std::endl;
}


void
test_parallel_fullydistributed()
{
  constexpr const unsigned dim = 2;

  parallel::distributed::Triangulation<dim> tria_base(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria_base);

  const auto description =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      tria_base, MPI_COMM_WORLD);

  parallel::fullydistributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  tria.create_triangulation(description);

  DoFHandler<dim> dh(tria);
  dh.distribute_dofs(FE_Q<dim>(1));

  tria.clear();
  deallog << "ndofs:" << dh.n_dofs() << std::endl;
}


int
main(int argc, char *argv[])
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  test_serial();
  test_parallel_shared();
  test_parallel_distributed();
  test_parallel_fullydistributed();
}
