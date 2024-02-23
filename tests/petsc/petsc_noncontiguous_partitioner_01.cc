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


// Test PETScWrappers::CommunicationPattern for non-contiguous index space.
// Copy-pasted from tests/base/mpi_noncontiguous_partitioner_01.cc

#include <deal.II/base/mpi.h>

#include <deal.II/lac/petsc_communication_pattern.h>

#include "../tests.h"


void
test(const MPI_Comm comm)
{
  IndexSet index_set_has(4);
  IndexSet index_set_want(4);

  if (Utilities::MPI::this_mpi_process(comm) == 0)
    {
      index_set_has.add_index(1);
      index_set_want.add_index(2);
    }
  else
    {
      index_set_has.add_index(2);
      index_set_want.add_index(1);
      index_set_want.add_index(2);
    }

  PETScWrappers::CommunicationPattern petscsf;
  petscsf.reinit(index_set_has, index_set_want, comm);

  AlignedVector<double> src(index_set_has.n_elements());
  AlignedVector<double> dst(index_set_want.n_elements());

  src[0] = Utilities::MPI::this_mpi_process(comm) * 100 + 1;

  petscsf.export_to_ghosted_array(ArrayView<const double>(src.data(),
                                                          src.size()),
                                  ArrayView<double>(dst.data(), dst.size()));

  for (size_t i = 0; i < src.size(); ++i)
    deallog << static_cast<int>(src[i]) << ' ';
  deallog << std::endl;
  for (size_t i = 0; i < dst.size(); ++i)
    deallog << static_cast<int>(dst[i]) << ' ';
  deallog << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  const MPI_Comm comm = MPI_COMM_WORLD;

  {
    deallog.push("all");
    test(comm);
    deallog.pop();
  }
}
