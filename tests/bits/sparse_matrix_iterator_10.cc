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



// this test, extracted from dof_constraints_09 and sparse_matrix_iterator_10,
// used to fail with aborts. this test is equivalent to
// sparse_matrix_iterator_09, except that we use the postfix operator++
// instead of the prefix form

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

  // and loop over the elements of it
  for (SparseMatrix<double>::const_iterator k = A.begin(); k != A.end(); ++k)
    deallog << k->row() << ' ' << k->column() << ' ' << k->value() << std::endl;
}



int
main()
{
  initlog();

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
