// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2017 by the deal.II authors
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



// like sparse_matrix_iterator_11, but for const_iterators

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

#include "../tests.h"


void
test()
{
  // create a sparsity pattern with totally
  // empty lines (not even diagonals, since
  // not quadratic)
  SparsityPattern sparsity(4, 5, 1);
  sparsity.add(1, 1);
  sparsity.add(3, 1);
  sparsity.compress();

  // attach a sparse matrix to it
  SparseMatrix<double> A(sparsity);

  SparseMatrix<double>::const_iterator k = A.begin(), j = ++A.begin();

  DEAL_II_AssertThrow(k < j, ExcInternalError());
  DEAL_II_AssertThrow(j > k, ExcInternalError());

  DEAL_II_AssertThrow(!(j < k), ExcInternalError());
  DEAL_II_AssertThrow(!(k > j), ExcInternalError());

  DEAL_II_AssertThrow(k != j, ExcInternalError());
  DEAL_II_AssertThrow(!(k == j), ExcInternalError());

  DEAL_II_AssertThrow(k == k, ExcInternalError());
  DEAL_II_AssertThrow(!(k != k), ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  try
    {
      test();
    }
  catch (std::exception &exc)
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
