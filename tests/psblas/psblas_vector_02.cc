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
#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/types.h>

#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/psblas_vector.h>

#include <psb_c_dbase.h>

#include <iostream>

#include "../../tests/tests.h"

using namespace dealii;

// Test the behavior of PSBLAS vectors with and without ghost elements.
// In particular:
// - create a PSBLAS vector without ghost elements and set it to a known
// non-ghosted vector
// - test a few norms
// - test range-based for loops(i.e. begin() and end())

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
    psblas_vector(idx) += idx;
  psblas_vector.compress(VectorOperation::add);

  PSCToolkitWrappers::Vector test_ghosted;
  test_ghosted.reinit(locally_owned_dofs,
                      locally_relevant_dofs,
                      mpi_communicator);
  test_ghosted = psblas_vector;
  AssertThrow(*test_ghosted.begin() == *psblas_vector.begin(),
              ExcMessage("First entry does not match!"));

  types::global_dof_index idx = 0;
  for (const double val : test_ghosted)
    {
      if (locally_owned_dofs.is_element(idx))
        Assert(val - psblas_vector(idx) < 1e-15,
               ExcMessage("Entries do not match!"));
      ++idx;
    }

  // let's clear the vector, using operator= when the vector is empty
  test_ghosted.clear();
  AssertThrow(test_ghosted.size() == 0, ExcInternalError());
  test_ghosted = psblas_vector;

  AssertThrow(test_ghosted.l1_norm() - psblas_vector.l1_norm() < 1e-15,
              ExcMessage("Norms do not match!"));
  AssertThrow(test_ghosted.l2_norm() - psblas_vector.l2_norm() < 1e-15,
              ExcMessage("Norms do not match!"));
  AssertThrow(test_ghosted.linfty_norm() - psblas_vector.linfty_norm() < 1e-15,
              ExcMessage("Norms do not match!"));

  // Test size()
  psblas_vector.clear();
  AssertThrow(psblas_vector.size() == 0, ExcInternalError());
  deallog << "OK" << std::endl;


  return 0;
}
