// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check copy_to() and copy_from() between tensors and full matrices.

#include <deal.II/base/tensor.h>

#include <deal.II/lac/full_matrix.h>

#include "../tests.h"

int
main()
{
  initlog();

  const unsigned int matrix_dimension = 4;

  Tensor<2, matrix_dimension> myTensor;

  FullMatrix<double> myMatrix(matrix_dimension, matrix_dimension);

  for (unsigned int i = 0; i < matrix_dimension; i++)
    for (unsigned int j = 0; j < matrix_dimension; j++)
      {
        myMatrix(i, j) = i + j;
      }

  myMatrix.copy_to(myTensor);

  deallog << matrix_dimension << std::endl;
  deallog.get_file_stream() << myTensor << std::endl;

  myTensor *= 2;
  myMatrix.copy_from(myTensor);
  for (unsigned int i = 0; i < matrix_dimension; i++)
    for (unsigned int j = 0; j < matrix_dimension; j++)
      {
        Assert(myMatrix(i, j) == 2. * (i + j), ExcInternalError());
      }
}
