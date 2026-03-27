// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Check that LinearAlgebra::TpetraWrappers::Vector::reinit
// does not drop ghost entries

#include <deal.II/base/utilities.h>

#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


template <typename Number>
void
test()
{
  IndexSet     parallel_partitioner_owned(10);
  IndexSet     parallel_partitioner_relevant(10);
  unsigned int rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (rank == 0)
    {
      parallel_partitioner_owned.add_range(0, 5);
      parallel_partitioner_relevant.add_range(5, 7);
    }
  else
    {
      parallel_partitioner_owned.add_range(5, 10);
      parallel_partitioner_relevant.add_range(3, 5);
    }

  parallel_partitioner_owned.compress();
  parallel_partitioner_relevant.compress();

  LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace::Default> a(
    parallel_partitioner_owned, parallel_partitioner_relevant, MPI_COMM_WORLD);

  LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace::Default> b;

  b.reinit(a);

  std::stringstream vector_output;
  a.print(vector_output);
  deallog << "Vector a:" << vector_output.str() << std::endl;

  // reset stream
  vector_output.str("");

  b.print(vector_output);
  deallog << "Vector b:" << vector_output.str() << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  MPILogInitAll log;

  test<double>();

  deallog << "OK" << std::endl;

  return 0;
}
