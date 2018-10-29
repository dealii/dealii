// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
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



// test PETScWrappers::MatrixBase::const_iterator::operator=

#include <deal.II/lac/petsc_sparse_matrix.h>

#include <iostream>

#include "../tests.h"


void
test()
{
  const types::global_dof_index n_rows = 5;
  PETScWrappers::SparseMatrix   matrix(n_rows, n_rows, 2);
  matrix.set(0, 0, 1.0);
  matrix.set(0, n_rows - 1, -1.0);
  for (unsigned int row_n = 1; row_n < n_rows; ++row_n)
    {
      matrix.set(row_n, row_n, 1.0);
      matrix.set(row_n, row_n - 1, -1.0);
    }
  matrix.compress(VectorOperation::insert);

  for (PETScWrappers::SparseMatrix::const_iterator iterator = matrix.begin();
       iterator != matrix.end();
       ++iterator)
    {
      // This is what we want to test.
      PETScWrappers::SparseMatrix::const_iterator duplicate = iterator;
      deallog << "row: " << duplicate->row() << std::endl;
      deallog << "column: " << duplicate->column() << std::endl;
      deallog << "value: " << duplicate->value() << std::endl;
    }

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
  catch (std::exception &exc)
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
