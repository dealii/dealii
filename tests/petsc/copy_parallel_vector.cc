// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check that we can correctly copy a distributed PETSc parallel vector onto a
// single deal.II vector.

#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>

#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/vector.h>

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

  deallog << "rank is: " << rank << std::endl;
  owned_indices.print(deallog);

  PETScWrappers::MPI::Vector petsc_vector(owned_indices, MPI_COMM_WORLD);
  for (const auto dof : owned_indices)
    {
      petsc_vector(dof) = double(dof);
    }

  Vector<double> deal_vector(petsc_vector);
  deallog << rank << ": deal vector size: " << deal_vector.size() << std::endl;

  for (const auto value : deal_vector)
    {
      deallog << value << std::endl;
    }
}
