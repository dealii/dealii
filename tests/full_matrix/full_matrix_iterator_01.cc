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



// like sparse_matrix_iterator_12, but for FullMatrix

#include <deal.II/lac/full_matrix.h>

#include "../tests.h"


template <typename IteratorType>
void
test()
{
  FullMatrix<double> A(3, 3);

  // test prefix operator
  const IteratorType k = A.begin(), j = std::next(A.begin());

  AssertThrow(k < j, ExcInternalError());
  AssertThrow(j > k, ExcInternalError());

  AssertThrow(k <= j, ExcInternalError());
  AssertThrow(j >= k, ExcInternalError());

  AssertThrow(!(j < k), ExcInternalError());
  AssertThrow(!(k > j), ExcInternalError());

  AssertThrow(k != j, ExcInternalError());
  AssertThrow(!(k == j), ExcInternalError());

  AssertThrow(k == k, ExcInternalError());
  AssertThrow(!(k != k), ExcInternalError());

  // test postfix operator
  IteratorType l = A.begin();
  IteratorType m = l++;

  AssertThrow(m == k, ExcInternalError());
  AssertThrow(l > m, ExcInternalError());
  AssertThrow(l->column() == k->column() + 1, ExcInternalError());
}



int
main()
{
  initlog();

  try
    {
      test<FullMatrix<double>::iterator>();
      test<FullMatrix<double>::const_iterator>();
      deallog << "OK" << std::endl;
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
