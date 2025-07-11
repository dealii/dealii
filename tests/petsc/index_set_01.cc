// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Verify that IndexSet::make_petsc_is() works

#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>

#include <deal.II/lac/petsc_compatibility.h>

#include <iostream>

#include "../tests.h"

void
test(const IndexSet &index_set)
{
  IS is = index_set.make_petsc_is();

  PetscInt size = 0;
  auto     ierr = ISGetSize(is, &size);
  AssertThrow(ierr == PETSC_SUCCESS, ExcPETScError(ierr));
  deallog << "size = " << size << std::endl;

  PetscInt local_size = 0;
  ierr                = ISGetLocalSize(is, &local_size);
  AssertThrow(ierr == PETSC_SUCCESS, ExcPETScError(ierr));
  deallog << "local size = " << local_size << std::endl;

  const PetscInt *local_entries = nullptr;
  ierr                          = ISGetIndices(is, &local_entries);
  AssertThrow(ierr == PETSC_SUCCESS, ExcPETScError(ierr));

  deallog << "entries =" << std::endl;
  for (PetscInt i = 0; i < local_size; ++i)
    deallog << "  " << local_entries[i] << std::endl;

  ierr = ISRestoreIndices(is, &local_entries);
  AssertThrow(ierr == PETSC_SUCCESS, ExcPETScError(ierr));
}


int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  deallog << "testing empty IndexSet" << std::endl;
  {
    IndexSet index_set(100);
    test(index_set);
  }
  deallog << std::endl;

  deallog << "testing an IndexSet with holes" << std::endl;
  {
    IndexSet index_set(100);
    index_set.add_range(10, 20);
    index_set.add_range(90, 100);
    test(index_set);
  }
  deallog << std::endl;

  deallog << "testing a complete IndexSet" << std::endl;
  {
    IndexSet index_set = complete_index_set(100);
    test(index_set);
  }
}
