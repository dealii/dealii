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



// check FullMatrix::matrix_scalar_product

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <vector>

#include "../tests.h"


void
test(Vector<double> &v, Vector<double> &w)
{
  FullMatrix<double> m(v.size(), v.size());
  for (unsigned int i = 0; i < m.m(); ++i)
    for (unsigned int j = 0; j < m.m(); ++j)
      m(i, j) = (i + 2 * j);

  for (unsigned int i = 0; i < v.size(); ++i)
    {
      v(i) = i;
      w(i) = i + 1;
    }

  v.compress();
  w.compress();

  // <w,Mv>
  const double s = m.matrix_scalar_product(w, v);

  // make sure we get the expected result
  for (unsigned int i = 0; i < v.size(); ++i)
    {
      AssertThrow(v(i) == i, ExcInternalError());
      AssertThrow(w(i) == i + 1, ExcInternalError());
    }

  double result = 0;
  for (unsigned int i = 0; i < m.m(); ++i)
    for (unsigned int j = 0; j < m.m(); ++j)
      result += (i + 2 * j) * j * (i + 1);

  AssertThrow(s == result, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  try
    {
      Vector<double> v(100);
      Vector<double> w(100);
      test(v, w);
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
