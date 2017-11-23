// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2017 by the deal.II Authors
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

// test DoFRenumbering::subdomain_wise with parallel::shared::Triangulation

#include "../tests.h"

// all include files you need here
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>

template<int dim>
void test()
{
  parallel::shared::Triangulation<dim>
  tria(MPI_COMM_WORLD,
       ::Triangulation<dim>::none,
       false,
       parallel::shared::Triangulation<dim>::partition_zorder);

  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);

  FE_Q<dim> fe(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  DoFRenumbering::subdomain_wise(dof_handler);

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "OK for " << dim << "d" << std::endl;
}

int main (int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
  MPILogInitAll log;

  test<2>();
  test<3>();

  return 0;
}
