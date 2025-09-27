// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check that we can correctly serialize a vector via boost.

#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>

#include <deal.II/lac/petsc_vector.h>

#include "../tests.h"


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_context(argc, argv, 1);
  MPILogInitAll                    mpi_log;

  const auto rank            = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const auto n_mpi_processes = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  constexpr types::global_dof_index dofs_per_process = 5;
  IndexSet owned_indices(dofs_per_process * n_mpi_processes);
  owned_indices.add_range(dofs_per_process * rank,
                          dofs_per_process * (rank + 1));
  owned_indices.compress();

  owned_indices.print(deallog);

  PETScWrappers::MPI::Vector petsc_vector(owned_indices, MPI_COMM_WORLD);
  for (const auto dof : owned_indices)
    petsc_vector(dof) = double(dof);

  petsc_vector.compress(VectorOperation::insert);

  std::stringstream               stream;
  boost::archive::binary_oarchive out_archive(stream);
  out_archive << petsc_vector;

  PETScWrappers::MPI::Vector petsc_copy(owned_indices, MPI_COMM_WORLD);

  boost::archive::binary_iarchive in_archive(stream);
  in_archive >> petsc_copy;

  petsc_copy.print(deallog.get_file_stream());

  AssertThrow(petsc_copy == petsc_vector, ExcInternalError());

  deallog << "OK" << std::endl;
}
