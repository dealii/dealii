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



// check PETScWrappers::SparseMatrix::frobenius_norm

#include <deal.II/lac/petsc_sparse_matrix.h>

#include <iostream>

#include "../tests.h"


void
test(PETScWrappers::SparseMatrix &m)
{
  // first set a few entries. count how many entries we have
  PetscScalar norm = 0;
  for (unsigned int i = 0; i < m.m(); ++i)
    for (unsigned int j = 0; j < m.m(); ++j)
      if ((i + 2 * j + 1) % 3 == 0)
        {
          m.set(i, j, std::complex<double>(0., i * j * .5 + .5));
          norm += (i * j * .5 + .5) * (i * j * .5 + .5);
        }
  norm = std::sqrt(norm);

  m.compress(VectorOperation::insert);

  // compare against the exact value of the l2-norm (max row-sum)
  deallog << m.frobenius_norm() << std::endl;
  AssertThrow(m.frobenius_norm() == norm, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();
  deallog.depth_console(0);

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        PETScWrappers::SparseMatrix m(5, 5, 3);
        test(m);
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
