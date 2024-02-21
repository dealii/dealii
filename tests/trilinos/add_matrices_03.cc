// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test SparseMatrix::add(factor, SparseMatrix) based on matrices of the same
// sparsity pattern

#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <iostream>

#include "../tests.h"


void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "numproc=" << numproc << std::endl;


  // each processor owns 3 indices
  IndexSet local_owned(numproc * 3);
  local_owned.add_range(myid * 3, myid * 3 + 3);

  // Create sparsity patterns
  TrilinosWrappers::SparsityPattern sp(local_owned, MPI_COMM_WORLD);

  for (unsigned int i = myid * 3; i < myid * 3 + 3; ++i)
    for (unsigned int j = 0; j < local_owned.size(); ++j)
      if ((i + j) % 2 == 1)
        {
          sp.add(i, j);
        }

  sp.compress();

  // create matrices by adding some elements into the respective positions
  TrilinosWrappers::SparseMatrix m1(sp), m2(sp);
  for (unsigned int i = myid * 3; i < myid * 3 + 3; ++i)
    for (unsigned int j = 0; j < local_owned.size(); ++j)
      if ((i + j) % 2 == 1)
        {
          m1.add(i, j, i + j);
          if (j % 2 == 0)
            m2.add(i, j, i + 2 * j + 1);
        }
  m1.compress(VectorOperation::add);
  m2.compress(VectorOperation::add);

  m1.add(2, m2);

  // Check for correctness of entries (all floating point comparisons should
  // be exact)
  for (unsigned int i = myid * 3; i < myid * 3 + 3; ++i)
    for (unsigned int j = 0; j < local_owned.size(); ++j)
      if ((i + j) % 2 == 1 && j % 2 == 0)
        {
          Assert(m1.el(i, j) == (double)i + j + 2 * i + 4 * j + 2,
                 ExcInternalError());
        }
      else if ((i + j) % 2 == 1)
        {
          Assert(m1.el(i, j) == (double)i + j, ExcInternalError());
        }
      else
        {
          Assert(m1.el(i, j) == 0., ExcInternalError());
        }

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  try
    {
      test();
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
