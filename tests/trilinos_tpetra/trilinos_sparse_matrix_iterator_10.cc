// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2024 by the deal.II authors
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



// this test, extracted from dof_constraints_09 and
// trilinos_sparse_matrix_iterator_10, used to fail with aborts. this test is
// equivalent to trilinos_sparse_matrix_iterator_09, except that we use the
// postfix operator++ instead of the prefix form

#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>
#include <deal.II/lac/trilinos_tpetra_sparsity_pattern.h>

#include "../tests.h"


void
test()
{
  // create a sparsity pattern with totally
  // empty lines (not even diagonals, since
  // not quadratic)
  LinearAlgebra::TpetraWrappers::SparsityPattern<MemorySpace::Default> sparsity(
    4, 5, 1);
  sparsity.add(1, 1);
  sparsity.add(3, 1);
  sparsity.compress();

  // attach a sparse matrix to it
  LinearAlgebra::TpetraWrappers::SparseMatrix<double, MemorySpace::Default> A(
    sparsity);

  // and loop over the elements of it
  for (LinearAlgebra::TpetraWrappers::
         SparseMatrix<double, MemorySpace::Default>::const_iterator k =
           A.begin();
       k != A.end();
       k++)
    deallog << k->row() << ' ' << k->column() << ' ' << k->value() << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  try
    {
      test();
    }
  catch (const std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
