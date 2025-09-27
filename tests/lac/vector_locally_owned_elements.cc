// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check that LinearAlgebra::TpetraWrappers::Vector::locally_owned_elements
// returns the correct index set.

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/trilinos_tpetra_vector.templates.h>

#include "../tests.h"

template <typename VectorType>
void
do_test()
{
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  parallel::distributed::Triangulation<2> triangulation(mpi_communicator);
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(4);

  const FE_Q<2> fe(1);
  DoFHandler<2> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
  IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  VectorType tpetra_vector;
  tpetra_vector.reinit(locally_owned_dofs,
                       locally_relevant_dofs,
                       mpi_communicator);

  Assert(tpetra_vector.locally_owned_elements() == locally_owned_dofs,
         ExcInternalError());

  deallog << "OK" << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  initlog();

  do_test<LinearAlgebra::TpetraWrappers::Vector<double>>();
}
