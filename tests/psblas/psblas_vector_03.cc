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

// Test the behavior of extract_subvector_to for PSBLAS vectors with and without
// ghost elements

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
  test_ghosted = psblas_vector;

  // Testing extract_subvector_to()

  // First, check we can extract locally owned indices
  const unsigned int                   n_indices_to_extract = 5;
  std::vector<types::global_dof_index> indices;
  int                                  i = 0;
  for (types::global_dof_index idx : locally_owned_dofs)
    if (i++ < n_indices_to_extract)
      indices.push_back(idx);

  std::vector<double> values;
  values.resize(indices.size());
  psblas_vector.extract_subvector_to(indices.begin(),
                                     indices.end(),
                                     values.begin());

  deallog << "Extracted values for owned indices: " << std::endl;
  for (const double val : values)
    deallog << val << std::endl;

  // ... now add some ghosted indices too to the set
  if (id == 0)
    {
      indices.push_back(15);
      indices.push_back(16);
    }
  else
    {
      indices.push_back(13);
      indices.push_back(14);
    }

  values.resize(indices.size());
  test_ghosted.extract_subvector_to(indices.begin(),
                                    indices.end(),
                                    values.begin());

  deallog << "Extracted values for owned + ghosted indices: " << std::endl;
  for (const double val : values)
    deallog << val << std::endl;

  return 0;
}
