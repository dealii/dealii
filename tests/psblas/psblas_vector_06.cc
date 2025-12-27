// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/logstream.h>
#include <deal.II/base/types.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/psblas_vector.h>

#include <psb_c_dbase.h>

#include <iostream>

#include "../tests.h"

using namespace dealii;

// Test VectorOperations::insert for PSBLAS vectors.

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPI_Comm                         mpi_communicator = MPI_COMM_WORLD;

  AssertThrow(Utilities::MPI::n_mpi_processes(mpi_communicator) == 2,
              ExcMessage("This test needs to be run with 2 MPI processes."));

  initlog();

  int id;
  MPI_Comm_rank(mpi_communicator, &id);
  IndexSet locally_owned_dofs(25);
  if (id == 0)
    locally_owned_dofs.add_range(0, 15);
  else if (id == 1)
    locally_owned_dofs.add_range(15, 25);


  IndexSet locally_relevant_dofs(25);
  locally_relevant_dofs = locally_owned_dofs;
  if (id == 0)
    locally_relevant_dofs.add_range(15, 17);
  else if (id == 1)
    locally_relevant_dofs.add_range(12, 15);

  PSCToolkitWrappers::Vector psblas_vector(locally_owned_dofs,
                                           mpi_communicator);

  for (const types::global_dof_index idx : locally_owned_dofs)
    psblas_vector(idx) = 1.;
  psblas_vector.compress(VectorOperation::insert);

  AssertThrow(std::all_of(psblas_vector.begin(),
                          psblas_vector.end(),
                          [](const double val) { return val == 1.0; }),
              ExcMessage(
                "Not all values are equal to 1 after insert operation."));

  deallog << "OK" << std::endl;

  return 0;
}
