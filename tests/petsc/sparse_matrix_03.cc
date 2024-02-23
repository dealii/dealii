// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check SparseMatrix::transpose

#include <deal.II/lac/petsc_sparse_matrix.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test()
{
  const unsigned int          s = 3;
  PETScWrappers::SparseMatrix m(s, s, s);

  const auto print = [&]() {
    for (unsigned int i = 0; i < m.m(); ++i)
      {
        for (unsigned int j = 0; j < m.n(); ++j)
          deallog << m(i, j) << ' ';
        deallog << std::endl;
      }
  };

  for (unsigned int i = 0; i < m.m(); ++i)
    for (unsigned int j = 0; j < m.n(); ++j)
      m.set(i, j, i + m.m() * j);

  m.compress(VectorOperation::insert);

  deallog << "before: " << std::endl;
  print();

  m.transpose();

  deallog << "after: " << std::endl;
  print();

  m.transpose();

  deallog << "back to original: " << std::endl;
  print();

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  test();
}
