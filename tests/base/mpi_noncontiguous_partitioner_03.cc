// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test Utilities::MPI::NoncontiguousPartitioner for padding.

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_noncontiguous_partitioner.h>

#include "../tests.h"


void
test(const MPI_Comm                       comm,
     std::vector<types::global_dof_index> index_set_has,
     std::vector<types::global_dof_index> index_set_want)
{
  Utilities::MPI::NoncontiguousPartitioner vector;
  vector.reinit(index_set_has, index_set_want, comm);

  AlignedVector<double> src(index_set_has.size(), 0);
  AlignedVector<double> dst(index_set_want.size(), 0);

  for (unsigned int i = 0; i < index_set_has.size(); ++i)
    src[i] = Utilities::MPI::this_mpi_process(comm) * 100 + i;

  vector.export_to_ghosted_array(ArrayView<const double>(src.data(),
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

  const unsigned int rank = Utilities::MPI::this_mpi_process(comm);

  {
    deallog.push("padding-non");

    if (rank == 0)
      test(comm, {0, 1, 2, 3}, {4, 5, 6, 7});
    else
      test(comm, {4, 5, 6, 7}, {0, 1, 2, 3});
    deallog.pop();
  }

  {
    deallog.push("padding-src");

    if (rank == 0)
      test(comm, {0, 1, numbers::invalid_dof_index, 2, 3}, {4, 5, 6, 7});
    else
      test(comm, {4, 5, 6, 7}, {0, 1, 2, 3});
    deallog.pop();
  }

  {
    deallog.push("padding-dst");

    if (rank == 0)
      test(comm, {0, 1, 2, 3}, {4, 5, numbers::invalid_dof_index, 6, 7});
    else
      test(comm, {4, 5, 6, 7}, {0, 1, 2, 3});
    deallog.pop();
  }

  {
    deallog.push("duplicates");

    if (rank == 0)
      test(comm, {0, 1, 2, 3}, {4, 4, 5, 6, 7, 6});
    else
      test(comm, {4, 5, 6, 7}, {0, 1, 2, 1, 1, 3});
    deallog.pop();
  }
}
