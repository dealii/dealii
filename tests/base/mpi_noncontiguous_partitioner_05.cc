// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test Utilities::MPI::NoncontiguousPartitioner for non-contiguous index space
// and multiple components.

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_noncontiguous_partitioner.h>
#include <deal.II/base/mpi_noncontiguous_partitioner.templates.h>

#include "../tests.h"


void
test(const MPI_Comm comm)
{
  const unsigned int n_components = 2;

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

  Utilities::MPI::NoncontiguousPartitioner vector(index_set_has,
                                                  index_set_want,
                                                  comm);

  AlignedVector<double> src(n_components * index_set_has.n_elements());
  AlignedVector<double> dst(n_components * index_set_want.n_elements());

  src[0] = Utilities::MPI::this_mpi_process(comm) * 100 + 1;
  src[1] = Utilities::MPI::this_mpi_process(comm) * 100 + 2;

  vector.export_to_ghosted_array<double, n_components>(
    ArrayView<const double>(src.data(), src.size()),
    ArrayView<double>(dst.data(), dst.size()));

  for (size_t i = 0; i < src.size(); ++i)
    deallog << static_cast<int>(src[i]) << ' ';
  deallog << std::endl;
  for (size_t i = 0; i < dst.size(); ++i)
    deallog << static_cast<int>(dst[i]) << ' ';
  deallog << std::endl;

  dst.fill(0.0);
  vector.export_to_ghosted_array<double, 0>(
    ArrayView<const double>(src.data(), src.size()),
    ArrayView<double>(dst.data(), dst.size()),
    n_components);

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
