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



// hp::DoFHandler::distribute_dofs failed
// if a processor only owns artificial cells

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    typename Triangulation<dim>::MeshSmoothing(
      Triangulation<dim>::smoothing_on_refinement |
      Triangulation<dim>::smoothing_on_coarsening));

  GridGenerator::hyper_cube(tria, -1.0, 1.0);
  tria.refine_global(1);

  const unsigned int max_degree = (dim == 2) ? 4 : 8;

  DoFHandler<dim>       dh(tria);
  hp::FECollection<dim> fe_collection;
  for (unsigned int degree = 1; degree <= max_degree; ++degree)
    fe_collection.push_back(FE_Q<dim>(degree));

  dh.distribute_dofs(fe_collection); // <-- fails

  // make sure no processor is
  // hanging
  MPI_Barrier(MPI_COMM_WORLD);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);


  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();

      deallog.push("2d");
      test<2>();
      deallog.pop();
      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    {
      test<2>();
      test<3>();
    }
}
