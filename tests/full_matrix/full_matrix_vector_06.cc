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



// check FullMatrix::matrix_norm_square

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <vector>

#include "../tests.h"


void
test(Vector<double> &v)
{
  FullMatrix<double> m(v.size(), v.size());
  for (unsigned int i = 0; i < m.m(); ++i)
    for (unsigned int j = 0; j < m.m(); ++j)
      m(i, j) = (i + 2 * j);

  for (unsigned int i = 0; i < v.size(); ++i)
    v(i) = i;

  v.compress();

  // <w,Mv>
  const double s = m.matrix_norm_square(v);

  // make sure we get the expected result
  for (unsigned int i = 0; i < v.size(); ++i)
    AssertThrow(v(i) == i, ExcInternalError());

  double result = 0;
  for (unsigned int i = 0; i < m.m(); ++i)
    for (unsigned int j = 0; j < m.m(); ++j)
      result += (i + 2 * j) * j * i;

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
      test(v);
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
