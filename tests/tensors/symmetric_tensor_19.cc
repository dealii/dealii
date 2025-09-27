// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check inversion of rank-4 tensor

#include <deal.II/base/symmetric_tensor.h>

#include "../tests.h"



template <int dim>
void
check(const SymmetricTensor<4, dim> &A)
{
  const SymmetricTensor<4, dim> B = invert(A);

  // check left inverse
  deallog << "    checking left inverse" << std::endl;
  const SymmetricTensor<4, dim> T_left = B * A;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = 0; l < dim; ++l)
          {
            deallog << "      " << A[i][j][k][l] << ' ' << B[i][j][k][l] << ' '
                    << T_left[i][j][k][l] << std::endl;

            AssertThrow(std::fabs(T_left[i][j][k][l] -
                                  identity_tensor<dim>()[i][j][k][l]) < 1e-10,
                        ExcInternalError());
          }

  // check left inverse
  deallog << "    checking right inverse" << std::endl;
  const SymmetricTensor<4, dim> T_right = A * B;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = 0; l < dim; ++l)
          {
            deallog << "      " << A[i][j][k][l] << ' ' << B[i][j][k][l] << ' '
                    << T_right[i][j][k][l] << std::endl;

            AssertThrow(std::fabs(T_right[i][j][k][l] -
                                  identity_tensor<dim>()[i][j][k][l]) < 1e-10,
                        ExcInternalError());
          }
}



template <int dim>
void
test()
{
  deallog << "dim=" << dim << std::endl;

  deallog << "  unit tensor" << std::endl;
  check(identity_tensor<dim>());

  // do something with a more complicated
  // tensor. make sure it is not
  // rank-deficient, so choose elements by
  // hand
  deallog << "  complicated tensor" << std::endl;
  SymmetricTensor<4, dim> A;
  switch (dim)
    {
      case 1:
        A[0][0][0][0] = 2;
        break;

      case 2:
        A[0][0][0][0] = 2;
        A[0][0][1][1] = 4;
        A[0][0][0][1] = 8;
        A[1][1][0][0] = 4;
        A[1][1][1][1] = 6;
        A[1][1][0][1] = 10;
        A[0][1][0][0] = 6;
        A[0][1][1][1] = 10;
        A[0][1][0][1] = 16;
        break;

      case 3:
        // I'm too lazy to code something
        // up by hand here
        for (unsigned int i = 0; i < 3; ++i)
          for (unsigned int j = 0; j < 3; ++j)
            for (unsigned int k = 0; k < 3; ++k)
              for (unsigned int l = 0; l < 3; ++l)
                A[i][j][k][l] = random_value<double>();
        break;

      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
  check(A);
}



int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test<1>();
  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}
