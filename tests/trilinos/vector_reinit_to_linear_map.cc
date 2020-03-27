// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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

// run TrilinosWrappers::MPI::Vector::reinit() in MPI environment
// and with an IndexSet that is_ascending_and_one_to_one to check
// that make_trilinos_map will construct the object with a Map that
// is linear

#include <deal.II/base/index_set.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/trilinos_vector.h>

#include "../tests.h"


void
test()
{
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  parallel::distributed::Triangulation<2> tria(mpi_communicator);

  GridGenerator::hyper_cube(tria, -1, 0);
  tria.refine_global(2);

  const unsigned int poly_degree = 1;
  FE_Q<2>            fe(poly_degree);

  DoFHandler<2> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();

  TrilinosWrappers::MPI::Vector vector_linear;
  vector_linear.reinit(locally_owned_dofs, mpi_communicator);

  const int is_linear_index_set =
    locally_owned_dofs.is_ascending_and_one_to_one(mpi_communicator);

  const int is_linear_map = vector_linear.trilinos_vector().Map().LinearMap();

  if (is_linear_index_set == 1 && is_linear_map == 1)
    {
      deallog << "OK" << std::endl;
    }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  mpi_initlog();

  test();
}
