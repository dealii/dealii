// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test Utilities::MPI::NoncontiguousPartitioner::import_from_ghosted_array().

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_noncontiguous_partitioner.h>

#include "../tests.h"


void
test(const MPI_Comm comm)
{
  std::vector<types::global_dof_index> index_set_has;
  std::vector<types::global_dof_index> index_set_want;

  if (Utilities::MPI::this_mpi_process(comm) == 0)
    {
      index_set_has.push_back(0);
      index_set_has.push_back(1);
      index_set_has.push_back(2);

      index_set_want.push_back(0);
      index_set_want.push_back(1);
      index_set_want.push_back(2);
      index_set_want.push_back(3);
      index_set_want.push_back(5);
    }
  else
    {
      index_set_has.push_back(3);
      index_set_has.push_back(4);
      index_set_has.push_back(5);

      index_set_want.push_back(0);
      index_set_want.push_back(2);
      index_set_want.push_back(3);
      index_set_want.push_back(4);
      index_set_want.push_back(5);
    }

  Utilities::MPI::NoncontiguousPartitioner vector(index_set_has,
                                                  index_set_want,
                                                  comm);

  AlignedVector<double> src(index_set_want.size());
  AlignedVector<double> dst(index_set_has.size());

  for (unsigned int i = 0; i < src.size(); ++i)
    src[i] = i + Utilities::MPI::this_mpi_process(comm) * 100 + 1;

  vector.import_from_ghosted_array(VectorOperation::add,
                                   ArrayView<double>(src.data(), src.size()),
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
