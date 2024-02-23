// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test the various extract_subvector_to functions for
// parallel vectors and block vectors

#include <deal.II/base/index_set.h>

#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/la_parallel_block_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


template <typename VectorType>
void
set(VectorType &vector)
{
  for (unsigned int i = 0; i < vector.size(); ++i)
    if (vector.locally_owned_elements().is_element(i))
      vector(i) = i;
  vector.compress(VectorOperation::insert);
}


template <typename VectorType>
void
test(VectorType &vector)
{
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  // select every other element
  std::vector<typename VectorType::size_type> indices;
  for (unsigned int j = 0; j < vector.size() / 2; ++j)
    indices.push_back(2 * j);

  // do the extraction with the function that takes indices, then
  // assert correctness
  std::vector<typename VectorType::value_type> values1(indices.size());
  vector.extract_subvector_to(indices, values1);
  for (unsigned int j = 0; j < vector.size() / 2; ++j)
    AssertThrow(get_real_assert_zero_imag(values1[j]) == 2 * j,
                ExcInternalError());

  // do the same with the version of the function that takes iterators
  std::vector<typename VectorType::value_type> values2(indices.size());
  vector.extract_subvector_to(indices.begin(), indices.end(), values2.begin());
  for (unsigned int j = 0; j < vector.size() / 2; ++j)
    AssertThrow(get_real_assert_zero_imag(values2[j]) == 2 * j,
                ExcInternalError());

  // done
  if (myid == 0)
    deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  {
    IndexSet local(10);
    if (myid == 0)
      local.add_range(0, 5);
    if (myid == 1)
      local.add_range(5, 10);

    IndexSet dense_local(10);
    dense_local.add_range(0, 10);

    {
      deallog.push("deal.II");
      LinearAlgebra::distributed::Vector<double> w(local, MPI_COMM_WORLD);
      set(w);
      LinearAlgebra::distributed::Vector<double> v(local,
                                                   dense_local,
                                                   MPI_COMM_WORLD);
      v = w; // get copy of vector including ghost elements
      test(v);
      deallog.pop();
    }

    {
      deallog.push("PETSc");
      PETScWrappers::MPI::Vector w(local, MPI_COMM_WORLD);
      set(w);
      PETScWrappers::MPI::Vector v(local, dense_local, MPI_COMM_WORLD);
      v = w; // get copy of vector including ghost elements
      test(v);
      deallog.pop();
    }

    {
      deallog.push("Trilinos");
      TrilinosWrappers::MPI::Vector w(local, MPI_COMM_WORLD);
      set(w);
      TrilinosWrappers::MPI::Vector v(local, dense_local, MPI_COMM_WORLD);
      v = w; // get copy of vector including ghost elements
      test(v);
      deallog.pop();
    }


    std::vector<IndexSet> partitioning;
    {
      IndexSet block1(10);
      if (myid == 0)
        block1.add_range(0, 7);
      if (myid == 1)
        block1.add_range(7, 10);

      IndexSet block2(6);
      if (myid == 0)
        block2.add_range(0, 2);
      if (myid == 1)
        block2.add_range(2, 6);

      partitioning.push_back(block1);
      partitioning.push_back(block2);
    }

    std::vector<IndexSet> dense_partitioning;
    {
      IndexSet block1(10);
      block1.add_range(0, 10);

      IndexSet block2(6);
      block2.add_range(0, 6);

      dense_partitioning.push_back(block1);
      dense_partitioning.push_back(block2);
    }


    {
      deallog.push("deal.II");
      LinearAlgebra::distributed::BlockVector<double> w(partitioning,
                                                        MPI_COMM_WORLD);
      set(w);
      LinearAlgebra::distributed::BlockVector<double> v(partitioning,
                                                        dense_partitioning,
                                                        MPI_COMM_WORLD);
      v = w; // get copy of vector including ghost elements
      test(v);
      deallog.pop();
    }

    {
      deallog.push("PETSc");
      PETScWrappers::MPI::BlockVector w(partitioning, MPI_COMM_WORLD);
      set(w);
      PETScWrappers::MPI::BlockVector v(partitioning,
                                        dense_partitioning,
                                        MPI_COMM_WORLD);
      v = w; // get copy of vector including ghost elements
      test(v);
      deallog.pop();
    }

    {
      deallog.push("Trilinos");
      TrilinosWrappers::MPI::BlockVector w(partitioning, MPI_COMM_WORLD);
      set(w);
      TrilinosWrappers::MPI::BlockVector v(partitioning,
                                           dense_partitioning,
                                           MPI_COMM_WORLD);
      v = w; // get copy of vector including ghost elements
      test(v);
      deallog.pop();
    }
  }
}
