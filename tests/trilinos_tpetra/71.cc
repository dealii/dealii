// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2004 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// check LinearAlgebra::TpetraWrappers::Vector<double>::min() and max()

#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_tpetra_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


template <typename VectorType>
void
test()
{
  IndexSet parallel_partitioner_distributed(10);
  IndexSet parallel_partitioner_ghosted(10);

  unsigned int rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (rank == 0)
    {
      parallel_partitioner_distributed.add_range(0, 5);
      parallel_partitioner_ghosted.add_range(0, 7);
    }
  else
    {
      parallel_partitioner_distributed.add_range(5, 10);
      parallel_partitioner_ghosted.add_range(3, 10);
    }

  parallel_partitioner_distributed.compress();
  parallel_partitioner_ghosted.compress();

  VectorType vector_distributed(parallel_partitioner_distributed,
                                MPI_COMM_WORLD);
  VectorType vector_ghosted(parallel_partitioner_distributed,
                            parallel_partitioner_ghosted,
                            MPI_COMM_WORLD);

  // set some elements of the vector
  for (unsigned int i = 0; i < vector_distributed.size(); ++i)
    {
      vector_distributed(i) = i;
    }
  vector_distributed.compress(VectorOperation::insert);
  vector_ghosted = vector_distributed;


  // then check the min/max
  const double eps =
    typeid(typename VectorType::value_type) == typeid(double) ? 1e-14 : 1e-5;

  Assert(std::fabs(vector_ghosted.min()) < eps, ExcInternalError());
  Assert(std::fabs(vector_ghosted.max() - (vector_ghosted.size() - 1)) < eps,
         ExcInternalError());
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
#ifdef DEAL_II_TRILINOS_WITH_TPETRA_INST_DOUBLE
        test<LinearAlgebra::TpetraWrappers::Vector<double,
                                                   MemorySpace::Default>>();
#endif
#ifdef DEAL_II_TRILINOS_WITH_TPETRA_INST_FLOAT
        test<
          LinearAlgebra::TpetraWrappers::Vector<float, MemorySpace::Default>>();
#endif
      }
      deallog << "OK" << std::endl;
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
}
