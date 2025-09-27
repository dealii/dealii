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


// check SparseMatrix::add(factor, other)

#include <deal.II/lac/petsc_sparse_matrix.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test()
{
  const unsigned int          s = 10;
  PETScWrappers::SparseMatrix m(s, s, s);
  for (unsigned int k = 0; k < m.m(); ++k)
    for (unsigned int l = 0; l <= k; ++l)
      m.set(k, l, k + 2 * l);

  m.compress(VectorOperation::insert);

  PETScWrappers::SparseMatrix m2(s, s, s);
  m2.set(0, 1, 5.0);
  for (unsigned int k = 0; k < m2.m(); ++k)
    m2.set(k, k, PetscScalar(1.0 + k, -1.0 - k));
  m2.compress(VectorOperation::insert);

  // we now print the matrix m, it is all real. Then print after
  // adding m2 which is complex, and then subtract m2 again to get the
  // original matrix back.

  deallog << "before: " << m(0, 1) << std::endl;
  for (unsigned int k = 0; k < s; ++k)
    deallog << m(k, k) << ' ';
  deallog << std::endl;

  m.add(1.0, m2);

  deallog << "after: " << m(0, 1) << std::endl;
  for (unsigned int k = 0; k < s; ++k)
    deallog << m(k, k) << ' ';
  deallog << std::endl;

  m.add(-1.0, m2);

  deallog << "back to original: " << m(0, 1) << std::endl;
  for (unsigned int k = 0; k < s; ++k)
    deallog << m(k, k) << ' ';
  deallog << std::endl;

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
