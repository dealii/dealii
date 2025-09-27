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



// check SparseMatrix::Tmmult

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test()
{
  // A = [1, 4, 7; 0, 1, 2; 0, 0, 1]
  PETScWrappers::SparseMatrix A(3, 3, 3);
  A.set(0, 0, 1);
  A.set(0, 1, 4);
  A.set(0, 2, 7);
  A.set(1, 1, 1);
  A.set(1, 2, 2);
  A.set(2, 2, 1);
  A.compress(VectorOperation::insert);

  // B = [1, 2, 3; 0, -3, -6; 0, 0, 0]
  PETScWrappers::SparseMatrix B(3, 3, 3);
  B.set(0, 0, 1);
  B.set(0, 1, 2);
  B.set(0, 2, 3);
  B.set(1, 1, -3);
  B.set(1, 2, -6);
  B.compress(VectorOperation::insert);

  // v = [2, 2, 2]
  IndexSet indices(3);
  indices.add_range(0, 3);
  PETScWrappers::MPI::Vector v(indices, MPI_COMM_WORLD);
  for (unsigned int i = 0; i < v.size(); ++i)
    v(i) = 2;
  v.compress(VectorOperation::insert);

  PETScWrappers::SparseMatrix C(3, 3, 3);

  // C := A^T diag(v) B = [2, 4, 6; 8, 10, 12; 14, 16, 18]
  A.Tmmult(C, B, v);

  // make sure we get the expected result
  for (unsigned int i = 0; i < C.m(); ++i)
    for (unsigned int j = 0; j < C.n(); ++j)
      AssertThrow(C(i, j) == 6 * i + 2 * j + 2, ExcInternalError());

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
        test();
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
