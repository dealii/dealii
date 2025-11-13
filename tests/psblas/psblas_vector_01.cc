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

#include <deal.II/base/exception_macros.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/types.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/psblas_vector.h>

#include <psb_c_dbase.h>

#include <iostream>

#include "../tests.h"

using namespace dealii;

// Test reinit for non-ghosted and ghosted PSBLAS vectors.

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPI_Comm                         mpi_communicator = MPI_COMM_WORLD;

  AssertThrow(Utilities::MPI::n_mpi_processes(mpi_communicator) == 2,
              ExcMessage("This test needs to be run with 2 MPI processes."));

  MPILogInitAll log;
  int           id;
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
    psblas_vector(idx) += idx;
  psblas_vector.compress(VectorOperation::add);

  PSCToolkitWrappers::Vector test_ghosted;
  test_ghosted.reinit(locally_owned_dofs,
                      locally_relevant_dofs,
                      mpi_communicator);
  test_ghosted = psblas_vector; // lhs has ghost elements, rhs does not

  AssertThrow(test_ghosted.l2_norm() == psblas_vector.l2_norm(),
              ExcInternalError());
  deallog << "OK" << std::endl;


  // Now let's test the case where the left hand side has a different size (like
  // 0)
  PSCToolkitWrappers::Vector test_empty;
  test_empty = psblas_vector;
  AssertThrow(test_empty.size() == psblas_vector.size(), ExcInternalError());
  AssertThrow(test_empty.l2_norm() == psblas_vector.l2_norm(),
              ExcInternalError());
  deallog << "OK" << std::endl;

  return 0;
}
