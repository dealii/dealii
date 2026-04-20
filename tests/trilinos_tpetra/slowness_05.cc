// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2005 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// this is part of a whole suite of tests that checks the relative speed of
// using Trilinos as compared to the speed of our own
// library. the tests therefore may not all actually use Trilinos, but they are
// meant to compare it
//
// this test compares element-wise vector operations in a random order

#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>

#include <iostream>

#include "../tests.h"

template <typename VectorType>
void
set(const std::vector<unsigned int> &permutation, VectorType &vector)
{
  const unsigned int N = vector.size();

  for (unsigned int i = 0; i < N; ++i)
    vector(permutation[i]) = i;
}

template <typename VectorType>
void
add(const std::vector<unsigned int> &permutation, VectorType &vector)
{
  const unsigned int N = vector.size();

  for (unsigned int i = 0; i < N; ++i)
    vector(permutation[i]) += i;
}

template <typename VectorType>
double
read(const std::vector<unsigned int> &permutation, const VectorType &vector)
{
  double             sum = 0.0;
  const unsigned int N   = vector.size();

  for (unsigned int i = 0; i < N; ++i)
    sum += static_cast<const double>(vector(permutation[i]));

  return sum;
}

template <typename VectorType>
void
test()
{
  // Vector entries
  const unsigned int N = 40000;

  // first find a random permutation of the
  // indices
  std::vector<unsigned int> permutation(N);
  {
    std::vector<unsigned int> unused_indices(N);
    for (unsigned int i = 0; i < N; ++i)
      unused_indices[i] = i;

    for (unsigned int i = 0; i < N; ++i)
      {
        // pick a random element among the
        // unused indices
        const unsigned int k = Testing::rand() % (N - i);
        permutation[i]       = unused_indices[k];

        // then swap this used element to the
        // end where we won't consider it any
        // more
        std::swap(unused_indices[k], unused_indices[N - i - 1]);
      }
  }

  // build the vector
  IndexSet   indices = complete_index_set(N);
  VectorType vector(indices, MPI_COMM_WORLD);

  // element-wise write access
  set(permutation, vector);

  // element-wise addition
  add(permutation, vector);

  // element-wise read
  const double sum = read(permutation, vector);

  deallog << "Vector norm: " << vector.l2_norm() << ". Vector sum: " << sum
          << std::endl;

  deallog << "Ok." << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  try
    {
      {
        test<
          LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Host>>();
        test<LinearAlgebra::distributed::Vector<double, MemorySpace::Host>>();
      }
    }
  catch (const std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
  return 0;
}
