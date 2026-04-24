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


// Check that various LinearAlgebra::TpetraWrappers::Vector initialization
// and assignment routines correctly set state variables.

#include <deal.II/base/utilities.h>

#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"

template <typename VectorType>
void
test_for_equality(VectorType a, VectorType b)
{
  Assert(a.has_ghost_elements() == b.has_ghost_elements(),
         ExcMessage("Difference in ghosted state."));
  Assert(a.is_compressed() == b.is_compressed(),
         ExcMessage("Difference in compressed state."));
  Assert(a.size() == b.size(), ExcMessage("Difference in global size."));
  Assert(a.locally_owned_size() == b.locally_owned_size(),
         ExcMessage("Difference in locally owned size."));
  Assert(a.locally_owned_elements() == b.locally_owned_elements(),
         ExcMessage("Difference in owned elements."));

  const unsigned int my_rank =
    dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  for (const auto i : a.locally_owned_elements())
    Assert(a(i) == b(i),
           ExcMessage("Difference on rank " + std::to_string(my_rank) +
                      " in vector element index <" + std::to_string(i) +
                      ">. Values: <" + std::to_string(a(i)) + "> and <" +
                      std::to_string(b(i)) + ">."));
}

template <typename VectorType>
void
test_vector_operations(VectorType a)
{
  // Check copy constructor
  VectorType b(a);

  test_for_equality(a, b);

  // Check clear()
  VectorType c;
  b.clear();
  test_for_equality(b, c);

  // Check assignment
  b = a;
  test_for_equality(a, b);

  // Check swap function
  VectorType d;
  b.swap(c);

  test_for_equality(a, c);
  test_for_equality(b, d);
}

template <typename VectorType>
void
test()
{
  IndexSet     parallel_partitioner_owned(5);
  IndexSet     parallel_partitioner_relevant(5);
  unsigned int rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (rank == 0)
    {
      parallel_partitioner_owned.add_range(0, 2);
      parallel_partitioner_relevant.add_range(0, 3);
    }
  else
    {
      parallel_partitioner_owned.add_range(2, 5);
      parallel_partitioner_relevant.add_range(1, 5);
    }

  parallel_partitioner_owned.compress();
  parallel_partitioner_relevant.compress();

  VectorType a(parallel_partitioner_owned, MPI_COMM_WORLD);

  // test on fully-distributd empty vector
  test_vector_operations(a);

  for (const auto i : a.locally_owned_elements())
    a(i) = i;

  // test on fully-distributed filled vector
  test_vector_operations(a);

  if (rank == 0)
    {
      a(2) += 2;
    }
  else
    {
      a(1) += 1;
    }

  // test on fully-distributd vector after compressing
  a.compress(VectorOperation::add);
  test_vector_operations(a);

  // test on ghosted vector
  VectorType b(parallel_partitioner_owned,
               parallel_partitioner_relevant,
               MPI_COMM_WORLD);
  b = a;
  test_vector_operations(b);
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  MPILogInitAll log;

  test<typename LinearAlgebra::TpetraWrappers::Vector<double,
                                                      MemorySpace::Host>>();

  deallog << "OK" << std::endl;

  return 0;
}
