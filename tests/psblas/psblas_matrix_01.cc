// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2004 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// check setting elements in a PSBLAS matrix

#include "deal.II/base/exception_macros.h"
#include "deal.II/base/exceptions.h"
#include "deal.II/base/index_set.h"

#include <deal.II/lac/psblas_sparse_matrix.h>

#include <iostream>

#include "../tests.h"


int
main(int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      MPILogInitAll                    log;
      {
        const IndexSet                      index_set = complete_index_set(5);
        PSCToolkitWrappers::SparsityPattern sparsity_pattern(index_set,
                                                             MPI_COMM_WORLD);
        // sparsity pattern
        unsigned int n_expected_elements = 0;
        for (unsigned int i = 0; i < 5; ++i)
          for (unsigned int j = 0; j < 5; ++j)
            if ((i + 2 * j + 1) % 3 == 0)
              {
                sparsity_pattern.add(i, j);
                n_expected_elements++;
              }

        sparsity_pattern.compress();

        // now set a few entries
        PSCToolkitWrappers::SparseMatrix matrix(sparsity_pattern,
                                                MPI_COMM_WORLD);
        for (unsigned int i = 0; i < matrix.m(); ++i)
          for (unsigned int j = 0; j < matrix.m(); ++j)
            if ((i + 2 * j + 1) % 3 == 0)
              matrix.set(i, j, i * j * .5 + .5);

        matrix.compress();

        deallog << "Matrix dimensions: " << matrix.m() << " x " << matrix.n()
                << std::endl;

        AssertThrow(
          matrix.n_nonzero_elements() == n_expected_elements,
          ExcMessage(
            "The number of nonzero elements in the matrix does not match the "
            "expected value."));
        deallog << "Nonzero elements: " << matrix.n_nonzero_elements()
                << std::endl;

        // check that the new elements are the right ones.
        for (unsigned int i = 0; i < matrix.m(); ++i)
          for (unsigned int j = 0; j < matrix.n(); ++j)
            if ((i + 2 * j + 1) % 3 == 0)
              {
                const double expected_value = i * j * .5 + .5;
                AssertThrow(matrix.el(i, j) == expected_value,
                            ExcMessage("Wrong matrix entries. Test failed."));
              }

        deallog << "OK" << std::endl;
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
}
