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



// check PETScWrappers::MatrixBase::clear ()

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test(PETScWrappers::MatrixBase &m)
{
  AssertThrow(m.m() != 0, ExcInternalError());
  AssertThrow(m.n() != 0, ExcInternalError());

  m.clear();

  AssertThrow(m.m() == 0, ExcInternalError());
  AssertThrow(m.n() == 0, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        PETScWrappers::SparseMatrix v(100, 100, 5);
        test(v);
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
