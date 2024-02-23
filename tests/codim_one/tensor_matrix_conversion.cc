// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// a short (a few lines) description of what the program does

#include "../tests.h"

// all include files you need here

#include <deal.II/base/tensor.h>

#include <deal.II/lac/full_matrix.h>

template <typename number>
void
fill_matrix(FullMatrix<number> &A)
{
  for (unsigned int i = 0; i < A.m(); ++i)
    for (unsigned int j = 0; j < A.n(); ++j)
      A(i, j) = number(i * A.n() + j + 1);
}

template <typename number>
void
display_matrix(FullMatrix<number> M)
{
  deallog << M.m() << 'x' << M.n() << " matrix" << std::endl;
  for (unsigned int i = 0; i < M.m(); ++i)
    {
      for (unsigned int j = 0; j < M.n(); ++j)
        deallog << M(i, j) << ' ';
      deallog << std::endl;
    }
}

template <int b>
void
fill_tensor_2(Tensor<2, b> &T)
{
  for (unsigned int i = 0; i < b; ++i)
    for (unsigned int j = 0; j < b; ++j)
      T[i][j] = i * b + j + 1;
}


template <int b>
void
display_tensor_2(Tensor<2, b> &T)
{
  deallog << b << 'x' << b << " tensor" << std::endl;
  for (unsigned int i = 0; i < b; ++i)
    {
      for (unsigned int j = 0; j < b; ++j)
        deallog << T[i][j] << ' ';
      deallog << std::endl;
    }
}

int
main()
{
  initlog();

  FullMatrix<double> A1(10, 10);
  Tensor<2, 3>       T1;
  fill_tensor_2(T1);

  for (unsigned int n = 0; n < 3; ++n)
    for (unsigned int i = 0; i < 10 - n; ++i)
      for (unsigned int j = 0; j < 10 - n; ++j)
        {
          A1.copy_from(T1, 0, n, 0, n, i, j);
          display_matrix(A1);
          A1 = 0;
        }

  FullMatrix<double> A2(3, 3);
  fill_matrix(A2);
  Tensor<2, 3> T2;
  for (unsigned int n = 0; n < 3; ++n)
    for (unsigned int i = 0; i < 3 - n; ++i)
      for (unsigned int j = 0; j < 3 - n; ++j)
        {
          A2.copy_to(T2, 0, n, 0, n, i, j);
          display_tensor_2(T2);
          T2 = 0;
        }



  return 0;
}
