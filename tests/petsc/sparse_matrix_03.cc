// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



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
