// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2021 by the deal.II authors
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
