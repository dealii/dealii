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



// document that SparseMatrix::operator= is broken
// now forbidden

#include <deal.II/lac/petsc_sparse_matrix.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test()
{
  PETScWrappers::SparseMatrix m(10, 10, 10);
  for (unsigned int i = 0; i < m.m(); ++i)
    for (unsigned int j = 0; j < m.m(); ++j)
      m.set(i, j, i + 2 * j);


  m.compress(VectorOperation::insert);

  {
    PETScWrappers::SparseMatrix m2;

    Mat mm  = m;
    Mat m2m = m2;

    // deallog << mm << ' ' << m2m << std::endl;

    // forbidden:  m2 = m;
    {
      // forbidden: PETScWrappers::SparseMatrix m3=m2;

      Mat mm  = m;
      Mat m2m = m2;
      // deallog << mm << ' ' << m2m << std::endl;
    }
  }
}


int
main(int argc, char **argv)
{
  initlog();

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        test();

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
