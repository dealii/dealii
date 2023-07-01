// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2018 by the deal.II authors
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


// deal.II includes
#include <deal.II/lac/petsc_sparse_matrix.h>

#include <cassert>
#include <iostream>

#include "../tests.h"

// sparse matrix elements
void
test(PETScWrappers::SparseMatrix &m)
{
  deallog << "Check matrix access" << std::endl;

  // fill up a matrix with some numbers
  for (unsigned int k = 0; k < m.m(); ++k)
    for (unsigned int l = 0; l < m.n(); ++l)
      if (k > l)
        m.set(k, l, PetscScalar(k + l, -1. * (k + l)));

  m.compress(VectorOperation::insert);

  // check the matrix is correctly filled
  for (unsigned int k = 0; k < m.m(); ++k)
    for (unsigned int l = 0; l < m.n(); ++l)
      AssertThrow(m(k, l).real() == -1. * (m(k, l).imag()), ExcInternalError());

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
        PETScWrappers::SparseMatrix m(10, 10, 10);
        test(m);

        deallog << "matrix:" << std::endl;
        m.print(deallog.get_file_stream());
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
    }

  deallog.get_file_stream() << std::endl;

  return 0;
}
